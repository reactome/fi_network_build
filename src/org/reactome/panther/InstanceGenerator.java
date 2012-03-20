/*
 * Created on Mar 13, 2006
 *
 */
package org.reactome.panther;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.xpath.XPath;

/**
 * This interface is used to generate GKInstance for the corresponding instances in the Panther
 * Pathways.
 * @author guanming
 */
public class InstanceGenerator {
    
    private PersistenceAdaptor reactomeAdaptor;
    
    public InstanceGenerator() {
    }
    
    public void setReactomeAdaptor(PersistenceAdaptor adaptor) {
        this.reactomeAdaptor = adaptor;
    }
    
    public GKInstance getCompartment(String compartmentID,
                                     String positionToCompartment,
                                     Map<String, GKInstance> idToInstanceMap) throws Exception {
        GKInstance rtn = idToInstanceMap.get(compartmentID);
        if (rtn == null)
            return rtn;
        if (positionToCompartment == null)
            return rtn;
        if (positionToCompartment.equals(PantherConstants.INSIDE))
            return rtn;
        String compName = rtn.getDisplayName();
        compName = positionToCompartment + " of " + compName;
        rtn = idToInstanceMap.get(compName);
        if (rtn != null)
            return rtn;
        rtn = searchInstance(ReactomeJavaConstants.EntityCompartment, compName);
        if (rtn != null) {
            // Try to speed up the searching
            idToInstanceMap.put(compName, rtn);
            return rtn;
        }
        // Have to construct one
        rtn = createInstance(ReactomeJavaConstants.EntityCompartment);
        rtn.setDisplayName(compName);
        rtn.addAttributeValue(ReactomeJavaConstants.name, compName);
        idToInstanceMap.put(compName, rtn);
        return rtn;
    }
    
    public GKInstance getRegulationType(String type) throws Exception {
        GKInstance regulationType = searchInstance(ReactomeJavaConstants.RegulationType, type);
        if (regulationType != null)
            return regulationType;
        regulationType = createInstance(ReactomeJavaConstants.RegulationType);
        regulationType.addAttributeValue(ReactomeJavaConstants.name, type);
        regulationType.setDisplayName(type);
        return regulationType;
    }
    
    public GKInstance createCompartment(Element compartmentElm, 
            PantherNameHandler nameHandler) throws Exception {
        String id = compartmentElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
        if (id == null)
            id = compartmentElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
//      id should NOT be null
        if (id.equals(PantherConstants.DEFAULT_COMPARTMENT_NAME))
            return null; // Don't need to do anything for a default compartment
        String longName = getCompartmentLongName(compartmentElm);
        // no celldesigner name is given. No need to convert the compartment information.
        if (longName != null)
            longName = nameHandler.cleanupName(longName);
        if (longName == null)
            // Use the default name
            longName = id;
        GKInstance instance = searchInstance(ReactomeJavaConstants.EntityCompartment, longName);
        if (instance != null)
            return instance;
        String name = compartmentElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
        instance = createInstance(ReactomeJavaConstants.EntityCompartment);
//      Use longname as displayName. For level 1 model, name is just like id in level 2
//      model, it is too short.
        instance.setDisplayName(longName);
        instance.addAttributeValue(ReactomeJavaConstants.name, longName);
        // Only level 2 name is meaningful
        if ((PantherConstants.SBML_NS != PantherConstants.SBML_NS_LEVEL1) &&
             !longName.equals(name))
            instance.addAttributeValue(ReactomeJavaConstants.name,
                    name);
//      Keep this information in defintion for debugging purpose.
        instance.setAttributeValue(ReactomeJavaConstants.definition, id);
        return instance;
    }
    
    private GKInstance searchInstance(String clsName, String displayName) throws Exception {
        Collection list = reactomeAdaptor.fetchInstanceByAttribute(clsName,
                "_displayName",
                "=",
                displayName);
        if (list == null || list.size() == 0)
            return null;
        GKInstance instance = (GKInstance) list.iterator().next();
        return instance;
    }
    
    private GKInstance searchInstance(String clsName, String attName, Object value) throws Exception {
        Collection list = reactomeAdaptor.fetchInstanceByAttribute(clsName,
                                               attName,
                                               "=",
                                               value);
        if (list == null || list.size() == 0)
            return null;
        return (GKInstance) list.iterator().next();
    }
    
    /**
     * 
     * @param clsName - Class name of instance to be created
     * @return GKInstance - generated Instance
     * @throws Exception
     */
    public GKInstance createInstance(String clsName) throws Exception {
        if (reactomeAdaptor instanceof XMLFileAdaptor) {
            XMLFileAdaptor fileAdaptor = (XMLFileAdaptor) reactomeAdaptor;
            return fileAdaptor.createNewInstance(clsName);
        }
        SchemaClass cls = reactomeAdaptor.getSchema().getClassByName(clsName);
        GKInstance instance = new GKInstance(cls);
        instance.setDbAdaptor(reactomeAdaptor);
        return instance;
    }
    
    @SuppressWarnings("unchecked")
    public GKInstance cloneEntity(GKInstance entity) throws Exception {
        GKInstance clone = createInstance(entity.getSchemClass().getName());
        for (Iterator it = entity.getSchemClass().getAttributes().iterator(); it.hasNext();) {
            SchemaAttribute att = (SchemaAttribute) it.next();
            List values = entity.getAttributeValuesList(att);
            if (att.getName().equals(ReactomeJavaConstants.DB_ID) ||
                att.getName().equals(ReactomeJavaConstants._displayName))
                continue;
            if (values == null || values.size() == 0)
                continue;
            clone.setAttributeValueNoCheck(att, new ArrayList(values));
        }
        return clone;
    }
    
    public GKInstance createEntity(Protein protein, String pathwayName) throws Exception {
        GKInstance newInstance = createEntity(PantherConstants.PROTEIN_TYPE_NAME);
        StringBuilder defintion = new StringBuilder();
        defintion.append(protein.getProteinId()).append(", ");
        defintion.append(PantherConstants.PROTEIN_TYPE_NAME);
        defintion.append(", ").append(protein.getType());
        defintion.append(", ").append(pathwayName);
        newInstance.setAttributeValue(ReactomeJavaConstants.definition, defintion.toString());
        newInstance.addAttributeValue(ReactomeJavaConstants.name, protein.getName());
        newInstance.setDisplayName(protein.getName());
        protein.addEntity(newInstance);
        return newInstance;
    }
    
    public GKInstance createSummation(String text) throws Exception {
        if (text == null || text.length() == 0)
            return null;
        GKInstance summation = createInstance(ReactomeJavaConstants.Summation);
        summation.setAttributeValue(ReactomeJavaConstants.text,
                                    text);
        summation.setDisplayName(InstanceDisplayNameGenerator.generateDisplayName(summation));
        return summation;
    }
    
    public boolean isSimpleType(String clsType) {
        if (clsType.equals(PantherConstants.ION_TYPE_NAME) ||
            clsType.equals(PantherConstants.SIMPLE_MOLECULE_TYPE_NAME) ||
            clsType.equals(PantherConstants.UNKNOWN_TYPE_NAME) ||
            clsType.equals(PantherConstants.PHENOTYPE_TYPE_NAME))
            return true;
        return false;
    }
    
    public boolean isVertexType(String clsType) {
    	if (clsType.equals(ReactomeJavaConstants.EntityVertex) ||
    		clsType.equals(ReactomeJavaConstants.PathwayVertex) ||
    		clsType.equals(ReactomeJavaConstants.ReactionVertex))
    		return true;
    	return false;
    }
    
    public boolean intersect(Collection list1, Collection list2) {
        if (list1 == null || list2 == null)
            return false;
        for (Iterator it = list1.iterator(); it.hasNext();) {
            if (list2.contains(it.next()))
                return true;
        }
        return false;
    }
    
    public GKInstance searchSimpleEntity(String clsType, 
                                         List<String> names, 
                                         GKInstance compartment)
                                         throws Exception {
        String clsName = null;
        if (clsType.equals(PantherConstants.ION_TYPE_NAME) ||
            clsType.equals(PantherConstants.SIMPLE_MOLECULE_TYPE_NAME))
            clsName = ReactomeJavaConstants.SimpleEntity;
        else if (clsType.equals(PantherConstants.PHENOTYPE_TYPE_NAME) ||
                 clsType.equals(PantherConstants.UNKNOWN_TYPE_NAME))
            clsName = ReactomeJavaConstants.OtherEntity;
        if (clsName == null)
            return null;
        SchemaClass cls = reactomeAdaptor.getSchema().getClassByName(clsName);
        Collection instances = reactomeAdaptor.fetchInstancesByClass(cls);
        GKInstance instance = null;
        for (Iterator it = instances.iterator(); it.hasNext();) {
            instance = (GKInstance) it.next();
            GKInstance value = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.compartment);
            if (value != compartment)
                continue;
            List values = instance.getAttributeValuesList(ReactomeJavaConstants.name);
            if (intersect(values, names))
                return instance;
        }
        return null;
    }
    
    public GKInstance createEntity(String clsType) throws Exception {
        GKInstance rtn = null;
        if (clsType.equals(PantherConstants.PROTEIN_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        }
        else if (clsType.equals(PantherConstants.GENE_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        }
        else if (clsType.equals(PantherConstants.RNA_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        }
        else if (clsType.equals(PantherConstants.ANTI_SENSE_RNA_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        }
        else if (clsType.equals(PantherConstants.ION_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.SimpleEntity);
        }
        else if (clsType.equals(PantherConstants.SIMPLE_MOLECULE_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.SimpleEntity);
        }
        else if (clsType.equals(PantherConstants.UNKNOWN_TYPE_NAME)) {
            rtn = createInstance(ReactomeJavaConstants.OtherEntity);
        }
        else if (clsType.equals(PantherConstants.PHENOTYPE_TYPE_NAME) ||
                clsType.equals(PantherConstants.DEGRADED_TYPE_NAME)) {
            // Use OtherEntity for these two types to keep the original network structures.
            // The manual cleanup should make necessary changes.
            rtn = createInstance(ReactomeJavaConstants.OtherEntity);
            // rtn = createInstance(ReactomeJavaConstants.Pathway);
            // Add to top pathway.
            // topPathway.addAttributeValue(ReactomeJavaConstants.hasComponent,
            //                             rtn);
        }
        else if (clsType.equals(PantherConstants.COMPLEX_TYPE_NAME))
            rtn = createInstance(ReactomeJavaConstants.Complex);
        return rtn;
    }
    
    /**
     * creates a PathwayDiagram Item to save coordinates
     * @param clsType ReactomeJavaConstant types for Vertexes, Edges and PathwayDiagram
     * @return GKInstance - created Instance
     * @throws Exception
     */    
    public GKInstance createPathwayDiagramItem(String clsType) throws Exception {
    	GKInstance rtn = null;
    	if (clsType.equals(ReactomeJavaConstants.PathwayDiagram) ||
    			clsType.equals(ReactomeJavaConstants.ReactionVertex) ||
    			clsType.equals(ReactomeJavaConstants.PathwayVertex) ||
    			clsType.equals(ReactomeJavaConstants.EntityVertex) || 
    			clsType.equals(ReactomeJavaConstants.Edge)) {
    		rtn = createInstance(clsType);
    	}
    	return rtn;
    }
    
    protected String getCompartmentLongName(Element compartmentElm) throws JDOMException {
        String path = ".//celldesigner:name";
        Element nameElm = (Element) XPath.selectSingleNode(compartmentElm, path);
        if (nameElm != null)
            return nameElm.getTextTrim();
        return null;
    }
    
    public GKInstance createLiteratureReference(String referenceLine) throws Exception {
        GKInstance lrInstance = null;
        if (referenceLine.startsWith(PantherConstants.MEDLINE_LABEL) ||
            referenceLine.startsWith(PantherConstants.PMID_LABEL)) {
            String id = null;
            if (referenceLine.startsWith(PantherConstants.MEDLINE_LABEL))
                id = referenceLine.substring(PantherConstants.MEDLINE_LABEL.length() + 1).trim();
            else
                id = referenceLine.substring(PantherConstants.PMID_LABEL.length() + 1).trim();
            lrInstance = searchInstance(ReactomeJavaConstants.LiteratureReference,
                                        ReactomeJavaConstants.pubMedIdentifier,
                                        new Integer(id));
            if (lrInstance == null) {
                lrInstance = createInstance(ReactomeJavaConstants.LiteratureReference);
                lrInstance.setAttributeValue(ReactomeJavaConstants.pubMedIdentifier, new Integer(id));
                lrInstance.setDisplayName(id);
            }
        }
        else if (referenceLine.indexOf(PantherConstants.PUBMED_LINE_MARK) > 0) {
            // Get the pubmed ids
            int index = referenceLine.indexOf(PantherConstants.PUBMED_ID_MARK);
            index += PantherConstants.PUBMED_ID_MARK.length();
            String id = referenceLine.substring(index);
            lrInstance = searchInstance(ReactomeJavaConstants.LiteratureReference,
                                        ReactomeJavaConstants.pubMedIdentifier,
                                        new Integer(id));
            if (lrInstance == null) {
                lrInstance = createInstance(ReactomeJavaConstants.LiteratureReference);
                lrInstance.setAttributeValue(ReactomeJavaConstants.pubMedIdentifier,
                        new Integer(id));
                // It seems not necessary to pull out all information from PubMED.
                lrInstance.setDisplayName(id); 
            }
        }
        else { 
            int index = 0;
            if (referenceLine.startsWith(PantherConstants.WEB_SITE_LABEL)) {
                index = PantherConstants.WEB_SITE_LABEL.length();
            }
            else if (referenceLine.startsWith(PantherConstants.FREE_TEXT_LABEL)) {
                index = PantherConstants.FREE_TEXT_LABEL.length();
            }
            String url = referenceLine.substring(index + 1); // Add 1 to remove "=".
            // Use journal to hold the url 
            lrInstance = searchInstance(ReactomeJavaConstants.LiteratureReference,
                                        ReactomeJavaConstants.journal,
                                        url);
            if (lrInstance == null) {
                lrInstance = createInstance(ReactomeJavaConstants.LiteratureReference);
                lrInstance.setAttributeValue(ReactomeJavaConstants.journal,
                        url); 
                lrInstance.setDisplayName(url);
            }
        }
        return lrInstance;
    }
    
    public GKInstance getPantherDBInstance() throws Exception {
        GKInstance dbInstance = searchInstance(ReactomeJavaConstants.ReferenceDatabase,
                                               PantherConstants.PANTHER_DB_NAME);
        
        if (dbInstance != null)
            return dbInstance;
        dbInstance = createInstance(ReactomeJavaConstants.ReferenceDatabase);
        dbInstance.setDisplayName(PantherConstants.PANTHER_DB_NAME);
        dbInstance.addAttributeValue(ReactomeJavaConstants.name,
                                     PantherConstants.PANTHER_DB_NAME);
        dbInstance.addAttributeValue(ReactomeJavaConstants.url,
                                     PantherConstants.PANTHER_URL);
        return dbInstance;
    }
    
    public GKInstance getDBIdentifier(String accession) throws Exception {
        GKInstance dbIDInstance = searchInstance(ReactomeJavaConstants.DatabaseIdentifier,
                                                 ReactomeJavaConstants.identifier,
                                                 accession);
        if (dbIDInstance != null)
            return dbIDInstance;
        dbIDInstance = createInstance(ReactomeJavaConstants.DatabaseIdentifier);
        // Need a database instance
        GKInstance dbInstance = getPantherDBInstance();
        dbIDInstance.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                       dbInstance);
        dbIDInstance.setAttributeValue(ReactomeJavaConstants.identifier, accession);
        InstanceDisplayNameGenerator.setDisplayName(dbIDInstance);
        return dbIDInstance;
    }
    
    public GKInstance getModificationInstance(Element modificationElm, Protein protein) throws Exception {
        String residue = modificationElm.getAttributeValue(PantherConstants.RESIDUE_ATT_NAME);
        String state = modificationElm.getAttributeValue(PantherConstants.STATE_ATT_NAME);
        String position = protein.getPosition(residue);
        String displayName = state;
        if (position != null)
            displayName = state + " at " + position;
        GKInstance modifiedResidue = searchInstance(ReactomeJavaConstants.ModifiedResidue, 
                                                    displayName);
        if (modifiedResidue == null) {
            modifiedResidue = createInstance(ReactomeJavaConstants.ModifiedResidue);
            // Keep all information in the _displayName for the time being
            modifiedResidue.setDisplayName(displayName);
        }
        return modifiedResidue;
    }
    
}
