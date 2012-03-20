/*
 * Created on Oct 15, 2008
 *
 */
package org.reactome.hibernate;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.hibernate.Query;
import org.hibernate.Session;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.junit.Test;
import org.reactome.fi.ReactomeAnalyzerTopicHelper;
import org.reactome.funcInt.FIAnnotation;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.funcInt.ReactomeSourceType;

/**
 * This query class is used to check the source of a FI. The source of a FI
 * should be something like: 
 *  1) any types from the interaction type attribute in the Interaction class
 *  2) complex for FIs extracted from a complex
 *  3) types for FIs extracted from a reaction:
 *      a). Input
 *      b). Catalysis
 *      c). Activation
 *      d). Inhibition
 *  4). FIs predicted based on NBC
 * @author wgm
 *
 */
public class FISourceTypeReader extends HibernateFIPersistence {
    private MySQLAdaptor dba;
    // Used to process GKInstance
    private ReactomeAnalyzerTopicHelper topicHelper;
    // Load types
    private Map<String, FIAnnotation> nameToType;
    private Map<String, String> typeToReverseTpeMap;
    
    public FISourceTypeReader() {
        topicHelper = new ReactomeAnalyzerTopicHelper();
    }
    
    /**
     * Initialize a source dba based on some hard-coded information.
     * @return
     * @throws Exception
     */
    public MySQLAdaptor initSourceDBA() throws Exception {
//        dba = new MySQLAdaptor("localhost",
//                               "reactome_plus_i_v2",
//                               "root",
//                               "macmysql01",
//                               3306);
        dba = new MySQLAdaptor("localhost",
                               "reactome_28_plus_i",
                               "root",
                               "macmysql01",
                               3306);
        return dba;
    }
    
    public void setSourceDBA(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    /**
     * Query the source type for a FI described by two protein or gene names.
     * TODO: Basically the following query is not right. The output FIs in names have been normalized
     * based on some standard mapping. Please see method HibernateFIReader.generateFIFileInGeneInHibernate.
     * Note: the normalization in the above mentioned method has been removed since a bug in the mapping
     * file from NCBI.
     * @param name1
     * @param name2
     * @return
     * @throws Exception
     */
    public Map<String, Set<String>> queryTypes(String name1,
                                               String name2,
                                               Session session) throws Exception {
        List<Interaction> interactions = new HibernateFIReader().queryFIsBasedOnGeneNames(name1, 
                                                                                          name2, 
                                                                                          session);
        if (interactions == null || interactions.size() == 0) {
            System.err.println(name1 + ", " + name2 + " cannot find an interaction!");
            return null;
        }
        // Make sure the correct name types
        Map<String, Set<String>> pairToTypes = new HashMap<String, Set<String>>();
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            Interaction interaction = (Interaction) it.next();
            Set<String> types = queryType(interaction);
            // Use the upper case only to avoid more than two pairs between two proteins
            // (e.g.): GNAS and PIK3R1
            String protein1 = interaction.getFirstProtein().getLabel().toUpperCase();
            String protein2 = interaction.getSecondProtein().getLabel().toUpperCase();
            String key = protein1 + "\t" + protein2;
            Set<String> set = pairToTypes.get(key);
            if (set == null) {
                set = new HashSet<String>();
                pairToTypes.put(key, set);
            }
            set.addAll(types);
        }
        return convertPassiveTypes(pairToTypes);
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
    public FIAnnotation queryType(String name1,
                            String name2,
                            Session session) throws Exception {
        if (nameToType == null)
            initInteractionTypes();
        List<Interaction> interactions = new HibernateFIReader().queryFIsBasedOnGeneNames(name1, 
                                                                                          name2, 
                                                                                          session);
        if (interactions == null || interactions.size() == 0) {
            System.err.println(name1 + ", " + name2 + " cannot find an interaction!");
            return null;
        }
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
            index = type.getDirection().indexOf("-");
            if (type.getDirection().length() > 1) {
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
    
    /**
     * Remove types ending with "_by"
     * @param pairToTypes
     * @return
     */
    private Map<String, Set<String>> convertPassiveTypes(Map<String, Set<String>> pairToTypes) {
        Map<String, Set<String>> convertedPairToTypes = new HashMap<String, Set<String>>();
        for (String pair : pairToTypes.keySet()) {
            Set<String> types = pairToTypes.get(pair);
            int index = pair.indexOf("\t");
            String name1 = pair.substring(0, index);
            String name2 = pair.substring(index + 1);
            String reverseKey = name2 + "\t" + name1;
            for (String type : types) {
                if (type.endsWith("_by")) {
                    String newType = null;
                    if (type.equals("inhibited_by")) {
                        newType = "inhibit";
                    }
                    else {
                        index = type.indexOf("d_by");
                        newType = type.substring(0, index);
                    }
                    Set<String> converted = convertedPairToTypes.get(reverseKey);
                    if (converted == null) {
                        converted = new HashSet<String>();
                        convertedPairToTypes.put(reverseKey, converted);
                    }
                    converted.add(newType);
                }
                else {
                    Set<String> original = convertedPairToTypes.get(pair);
                    if (original == null) {
                        original = new HashSet<String>();
                        convertedPairToTypes.put(pair, original);
                    }
                    original.add(type);
                }
            }
        }
        return convertedPairToTypes;
    }
    
    /**
     * Query a type for a passed Interaction.
     * @param interaction
     * @return
     * @throws Exception
     */
    private Set<String> queryType(Interaction interaction) throws Exception {
        if (interaction == null)
            return null;
        Set<String> types = new HashSet<String>();
        if (interaction.getEvidence() != null) {
            types.add("predicted");
            return types; // Predicted based on NBC
        }
        // Need to query the Reactome Source
        Set<ReactomeSource> sources = interaction.getReactomeSources();
        if (dba == null)
            initSourceDBA();
        for (ReactomeSource src : sources) {
            if (src.getSourceType() == ReactomeSourceType.COMPLEX)
                types.add("complex");
            else if (src.getSourceType() == ReactomeSourceType.INTERACTION) {
                String type = extractInteractionTypeInString(src);
                types.add(type);
            }
            else if (src.getSourceType() == ReactomeSourceType.REACTION) {
                String type = extractTypeFromReaction(src, interaction);
                types.add(type);
            }
            else if (src.getSourceType() == ReactomeSourceType.TARGETED_INTERACTION) {
                String type = extractTypeFromTargetedInteraction(src, interaction);
                types.add(type);
            }
        }
        return types;
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
        if (dba == null)
            initSourceDBA();
        boolean needReverse = false;
        // Need to see if a reverse type should be used.
        if (interaction.getFirstProtein().getShortName().equals(name2))
            needReverse = true;
        FIAnnotation fiType = null;
        for (ReactomeSource src : sources) {
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
                if (needReverse)
                    fiType = fiType.generateReverseType();
            }
            else if (src.getSourceType() == ReactomeSourceType.TARGETED_INTERACTION) {
                String type = extractTypeFromTargetedInteraction(src, interaction);
                fiType = nameToType.get(type);
                if (needReverse)
                    fiType = fiType.generateReverseType();
            }
            if (fiType == null) {
                System.err.println(name1 + " " + name2 + " has no type!" + " Interaction: " + interaction.getDbId());
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
        if (src.getDataSource().equals("TRED")) { // Right now covers TRED TF/Target interactions only.
            Long dbId = src.getReactomeId();
            GKInstance instance = dba.fetchInstance(dbId);
            // Get factor
            GKInstance factor = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.factor);
            Set<GKInstance> refPepSeqs = topicHelper.grepRefPepSeqs(factor);
            Set<String> factorIds = new HashSet<String>();
            for (GKInstance refPepSeq : refPepSeqs) {
                String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                factorIds.add(id);
            }
            // Get target
            GKInstance target = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.target);
            refPepSeqs = topicHelper.grepRefPepSeqs(target);
            Set<String> targetIds = new HashSet<String>();
            for (GKInstance refPepSeq : refPepSeqs) {
                String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                targetIds.add(id);
            }
            String proteinId1 = interaction.getFirstProtein().getPrimaryAccession();
            String proteinId2 = interaction.getSecondProtein().getPrimaryAccession();
            if (factorIds.contains(proteinId1) && targetIds.contains(proteinId2)) {
                return "expression regulates";
            }
            if (factorIds.contains(proteinId2) && targetIds.contains(proteinId1)) {
                return "expression regulated by";
            }
        }
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
            Set<GKInstance> refPepSeqs = topicHelper.grepRefPepSeqs(instance);
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
        List interactors = instance.getAttributeValuesList(ReactomeJavaConstants.interactor);
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
     * Use the names to query proteins first. After proteins found, use found proteins
     * to query Interaction. Since we have not stored labels for proteins in the database,
     * we have to go through this detoured way.
     * @param name1
     * @param name2
     * @param session
     * @return
     * @throws Exception
     */
    private List queryFIBasedOnProteins(String name1,
                                        String name2,
                                        Session session) throws Exception {
        Query query = session.createQuery("FROM Protein as p where p.shortName = ? or p.name = ?");
        // Search the first protein
        query.setString(0, name1);
        query.setString(1, name1);
        List proteins1 = query.list();
        // Search the second protein
        query.setString(0, name2);
        query.setString(1, name2);
        List proteins2 = query.list();
        query = session.createQuery("FROM Interaction as i WHERE i.firstProtein = ? and i.secondProtein = ?");
        for (Iterator it1 = proteins1.iterator(); it1.hasNext();) {
            Protein protein1 = (Protein) it1.next();
            for (Iterator it2 = proteins2.iterator(); it2.hasNext();) {
                Protein protein2 = (Protein) it2.next();
                query.setEntity(0, protein1);
                query.setEntity(1, protein2);
                List interactions = query.list();
                if (interactions != null && interactions.size() > 0)
                    return interactions;
                // Try another way
                query.setEntity(0, protein2);
                query.setEntity(1, protein1);
                interactions = query.list();
                if (interactions != null && interactions.size() > 0)
                    return interactions;
            }
        }
        return null;
    }
    
    @Test
    public void testQueryProtein() throws Exception {
        initSession();
        Session session = sessionFactory.openSession();
        String name = "PAK3";
        Query query = session.createQuery("From Protein as p where p.shortName = ? or p.name = ?");
        query.setString(0, name);
        query.setString(1, name);
        List list = query.list();
        System.out.println("Total proteins: " + list.size());
    }
    
    /**
     * This method is used to list interaction types from the database.
     * @throws Exception
     */
    @Test
    public void listInteractionTypes() throws Exception {
        initSourceDBA();
        Collection interactions = dba.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Interaction);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.interactionType);
        dba.loadInstanceAttributeValues(interactions, att);
        Set<String> types = new HashSet<String>();
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            GKInstance interaction = (GKInstance) it.next();
            String type = (String) interaction.getAttributeValue(ReactomeJavaConstants.interactionType);
            if (type != null)
                types.add(type);
        }
        System.out.println("Types: " + types.size());
        for (String type : types) {// For the InteractionTypeMapper file.
            System.out.println("<type name=\"" + type + "\" reverse=\"" + type + "\" direction=\"none\"/>");
            //System.out.println(type);
        }
    }
    
    /**
     * Load FI types from a pre-generated file and add some types known already (e.g. predicted).
     * @return
     * @throws Exception
     */
    private void initInteractionTypes() throws Exception {
        String fileName = "resources/InteractionTypeMapper.xml";
        nameToType= new HashMap<String, FIAnnotation>();
        typeToReverseTpeMap = new HashMap<String, String>();
        SAXBuilder builder = new SAXBuilder();
        File file = new File(fileName);
        Document document = builder.build(file);
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
            typeToReverseTpeMap.put(name, reverse);
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
    
}
