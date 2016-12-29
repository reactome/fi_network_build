/*
 * Created on Feb 28, 2006
 *
 */
package org.reactome.panther;

import static org.reactome.panther.PantherConstants.MODEL_ELM_NAME;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.jaxen.JaxenException;
import org.jaxen.XPath;
import org.jdom.Attribute;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.reactome.convert.common.Converter;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.XPathHelper;

/**
 * This class is used to convert Panther pathways to Reactome pathways. Please note that only elements in
 * celldesigner namespace are processed in this class to take use of some properties in CellDesigner. So this
 * is not a generic converter from SBML2 to Reactome.<br><br>
 * 
 * This class will now import coordinate information from celldesigner. because of quite large differences
 * between some of the celldesigner-versions, this importer will only handle the versions 2.2 and 2.3. For the
 * import of versions 2.5 and 4.0, please go to <code>PantherToReactomeConverterV25</code>.<br>
 * 
 * @author guanming
 * @author andreash
 */
@SuppressWarnings("unchecked")
public class PantherToReactomeConverter implements Converter {
	static Logger logger = Logger.getLogger(PantherToReactomeConverter.class);
	// Used to load the Reactome Schema
	private XMLFileAdaptor fileAdaptor;
	// Database used to load housekeeping instances
	private MySQLAdaptor dbAdaptor;
	// State properties
	private Map<String, GKInstance> idToInstanceMap;
	// proteinReferences are used for complex definitions
	private Map<String, Protein> proteinMap;
	// Name are used for complex definitions for species that are not proteins.
	private Map<String, GKInstance> nameToInstanceMap;
	// Top level Pathway
	private GKInstance topPathway;
	// To avoid creating multiple CatalystActvities for the same physicalEntity
	private Map<GKInstance, GKInstance> caMap;
	// A helper class to handle names issues
	private PantherNameHandler nameHandler;
	// Used to control GKInstance generation
	private InstanceGenerator instanceGenerator;
	// degraded species uses meaningless name. Have to replace them
	private List<GKInstance> degradedEntities;
	// coordinates for aliases, global to avoid copying large data internally
	private Map<String, Map<String, String>> coordinates;
	// saves reaction and entity ids for later use in coordiante processing
	private Map<String, GKInstance> idToReactionMap;

	// maps reactionIDs (r1) to a map of functionality (input, output, catalyst)
	// and saves which one of the inputs/outputs are main members of the reaction
	// in mainInput and mainOutput
	private Map<String, Map<String, List<String>>> reactionToAliasMap;

	public PantherToReactomeConverter() {
		nameHandler = new PantherNameHandler();
		instanceGenerator = new InstanceGenerator();
	}

	public void setDatabaseAdaptor(MySQLAdaptor dbAdaptor) {
		this.dbAdaptor = dbAdaptor;
	}

	public MySQLAdaptor getDatabaseAdaptor() {
		return this.dbAdaptor;
	}

	public void setFileAdaptor(XMLFileAdaptor adaptor) {
		this.fileAdaptor = adaptor;
		instanceGenerator.setReactomeAdaptor(fileAdaptor);
	}

	public XMLFileAdaptor getFileAdaptor() {
		return this.fileAdaptor;
	}

	private void reset() {
		if (idToInstanceMap == null)
			idToInstanceMap = new HashMap<String, GKInstance>();
		else
			idToInstanceMap.clear();
		if (proteinMap == null)
			proteinMap = new HashMap<String, Protein>();
		else
			proteinMap.clear();
		if (nameToInstanceMap == null)
			nameToInstanceMap = new HashMap<String, GKInstance>();
		else
			nameToInstanceMap.clear();
		if (caMap == null)
			caMap = new HashMap<GKInstance, GKInstance>();
		else
			caMap.clear();
		if (degradedEntities == null)
			degradedEntities = new ArrayList<GKInstance>();
		else
			degradedEntities.clear();
		topPathway = null;

		if (idToReactionMap == null)
			idToReactionMap = new HashMap<String, GKInstance>();
		else
			idToReactionMap.clear();

		if (reactionToAliasMap == null)
			reactionToAliasMap = new HashMap<String, Map<String, List<String>>>();
		else
			reactionToAliasMap.clear();
	}

	public void setInstanceGenerator(InstanceGenerator instanceGenerator) {
		this.instanceGenerator = instanceGenerator;
	}

	public InstanceGenerator getInstanceGenerator() {
		return this.instanceGenerator;
	}

	/**
	 * Convert a list of SBML files into the pre-specified local project, which is specified by XMLFileAdaptor.
	 * @param fileNames
	 * @throws Exception
	 */
	public List<GKInstance> convert(List<String> fileNames) throws Exception {
		if (fileAdaptor == null)
			throw new IllegalStateException("PantherToReactomeConverter.convert(): No Reactome "
			      + "PersistenceAdaptor defined.");
		if (instanceGenerator == null)
			throw new IllegalStateException(
			      "PantherToReactomeConverer.convert(): No InstanceGenerator specified.");
		List<GKInstance> convertedPathways = new ArrayList<GKInstance>();
		for (String fileName : fileNames) {
			logger.info("Start converting " + fileName + "...");
			GKInstance pathway = convertBeforePost(fileName);
			convertedPathways.add(pathway);
		}
		logger.info("Start post-processing...");
		postProcess();
		return convertedPathways;
	}

	/**
	 * Convert a Panther SBML file to a local Reactome project.
	 * @param fileName SBML model Document
	 * @return the top GKInstance object for the converted model.
	 * @throws Exception
	 */
	public GKInstance convert(String fileName) throws Exception {
		if (fileAdaptor == null)
			throw new IllegalStateException(
			      "PantherToReactomeConverter.convert(): No Reactome PersistenceAdaptor defined.");
		if (instanceGenerator == null)
			throw new IllegalStateException(
			      "PantherToReactomeConverer.convert(): No InstanceGenerator specified.");
		GKInstance rPathway = convertBeforePost(fileName);
		postProcess();
		return rPathway;
	}

	public GKInstance convertBeforePost(String fileName) throws JDOMException, IOException, Exception {
		reset();
		// Load XML Document
		SAXBuilder builder = new SAXBuilder();
		Document sbmlDoc = builder.build(new File(fileName));
		Element root = sbmlDoc.getRootElement();
		// Some models are in level1. Level2 and Level2 have different
		// namespaces. Set SBML_NS dynamically.
		PantherConstants.SBML_NS = root.getNamespace();
		Element modelElm = root.getChild(MODEL_ELM_NAME, PantherConstants.SBML_NS);
		// Assume it works with XMLFileAdaptor
		GKInstance rPathway = createInstance(ReactomeJavaConstants.Pathway);
		extractPathwayProp(modelElm, rPathway);
		topPathway = rPathway;
		proteinMap = extractProteins(modelElm);
		// Get the SBML Model Version from the SBML File
		Element annotationElm = modelElm.getChild(PantherConstants.ANNOTATION_ELM_NAME,
		      PantherConstants.SBML_NS);
		// Handle compartments first since they are going to be used by species
		Element compartmentsElm = modelElm.getChild(PantherConstants.LIST_OF_COMPARTMENTS_ELM_NAME,
		      PantherConstants.SBML_NS);
		processCompartmentsElement(compartmentsElm);
		// Handle species
		Element speciesListElm = modelElm.getChild(PantherConstants.LIST_OF_SPECIES_ELM_NAME,
		      PantherConstants.SBML_NS);
		// SpeciesListElm should not be null
		processSpeciesElement(speciesListElm);

		// Some species have active states. PhysicalEntities should be created
		// for active states. Species aliases used in CellDesigner are mapped to
		// SBML species ids too.
		handleActiveStates(modelElm);
		// Check to make sure all complexes have names. Level1 pathways don't
		// have names
		nameHandler.validateComplexNames(idToInstanceMap);
		// Assign a meaningful name to degraded entities
		nameHandler.assignDegradedEntityNames(degradedEntities, idToInstanceMap);
		// Now it is the time for reactions
		Element reactionsElm = modelElm.getChild(PantherConstants.LIST_OF_REACTION_ELM_NAME,
		      PantherConstants.SBML_NS);
		processReactionsElement(reactionsElm);
		sanityCheckEvents();
		PostProcessHelper.generatePrecedingProperties(topPathway);

		//extractDiagram(rPathway, annotationElm, reactionsElm);
		return rPathway;
	}

    private void extractDiagram(GKInstance rPathway, 
                                Element annotationElm,
                                Element reactionsElm) throws  Exception {
        // handling of coordinates and vertices starts here
		extractCoordinates(annotationElm);
		GKInstance pathwayDiagram = createPathwayDiagram(rPathway);
		extractAliasMatching(reactionsElm);
		generateVertices(pathwayDiagram);
    }

	private void postProcess() throws Exception {
		if (dbAdaptor == null)
			throw new IllegalStateException(
			      "PantherToReactomeConverter.postProcess(): Database Adaptor is not specified!");
		PantherPostProcessor postProcessor = new PantherPostProcessor();
		postProcessor.postProcess(dbAdaptor, fileAdaptor);
		postProcessor.otherProcesses(dbAdaptor, fileAdaptor);
	}

	private List<GKInstance> searchInhibitedReactions(GKInstance input) throws Exception {
		List reactions = getReactions();
		List<GKInstance> rtn = new ArrayList<GKInstance>();
		GKInstance rxt = null;
		List inputs = null;
		List cas = null;
		for (Iterator it = reactions.iterator(); it.hasNext();) {
			rxt = (GKInstance) it.next();
			inputs = rxt.getAttributeValuesList(ReactomeJavaConstants.input);
			if (inputs != null && inputs.contains(input)) {
				rtn.add(rxt);
				continue;
			}
			cas = rxt.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
			if (cas != null && cas.size() > 0) {
				for (Iterator it1 = cas.iterator(); it1.hasNext();) {
					GKInstance ca = (GKInstance) it1.next();
					GKInstance catalyst = (GKInstance) ca.getAttributeValue(ReactomeJavaConstants.physicalEntity);
					if (input == catalyst) {
						rtn.add(rxt);
						break;
					}
				}
			}
		}
		return rtn;
	}

	/**
	 * Use this method so that this converter can work with both old and new schema: in the latest schema,
	 * hasComponent in pathway has been changed to hasEvent.
	 * @return
	 * @throws Exception
	 */
	private List getReactions() throws Exception {
		if (topPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
			return topPathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
		if (topPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent))
			return topPathway.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
		return null; // Should not be here!
	}

	private void addReaction(GKInstance reaction) throws Exception {
		if (topPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent)) {
			topPathway.addAttributeValue(ReactomeJavaConstants.hasEvent, reaction);
		} else if (topPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent)) {
			topPathway.addAttributeValue(ReactomeJavaConstants.hasComponent, reaction);
		}
	}

	/**
	 * A sanity check for events (reactions and pathways) is added: if a PhysicalEntity is in both input and
	 * output, take this PhysicalEntity out, create a Required instance for the Event instance.
	 * @throws Exception
	 */
	private void sanityCheckEvents() throws Exception {
		List events = getReactions();
		if (events == null || events.size() == 0)
			return;
		GKInstance event = null;
		for (Iterator it = events.iterator(); it.hasNext();) {
			event = (GKInstance) it.next();
			List inputs = event.getAttributeValuesList(ReactomeJavaConstants.input);
			List outputs = event.getAttributeValuesList(ReactomeJavaConstants.output);
			if (inputs == null || outputs == null)
				continue;
			for (Iterator it1 = inputs.iterator(); it1.hasNext();) {
				GKInstance input = (GKInstance) it1.next();
				if (outputs.contains(input)) {
					createRequirement(event, input);
					it1.remove();
					outputs.remove(input);
				}
			}
		}
	}

	private void createRequirement(GKInstance event, GKInstance regulator) throws Exception {
		GKInstance requirement = createInstance(ReactomeJavaConstants.Requirement);
		requirement.addAttributeValue(ReactomeJavaConstants.regulator, regulator);
		requirement.addAttributeValue(ReactomeJavaConstants.regulatedEntity, event);
		InstanceDisplayNameGenerator.setDisplayName(requirement);
	}

	/**
	 * Some species have active states. A new PhysicalEntity should be created for an active state.
	 * @param modelElm
	 * @throws Exception
	 */
	private void handleActiveStates(Element modelElm) throws Exception {
		// Get the compartment aliases first since it will be used for aliases
		Map<String, GKInstance> compartAliasMap = extractCompartmentAliases(modelElm);
		String xPath = ".//celldesigner:speciesAlias";
		XPath myPath = XPathHelper.getXPathObject(xPath, modelElm);
		List elements = myPath.selectNodes(modelElm);
		if (elements == null || elements.size() == 0)
			return;
		Map<GKInstance, GKInstance> activeMap = new HashMap<GKInstance, GKInstance>();
		for (Iterator it = elements.iterator(); it.hasNext();) {
			Element elm = (Element) it.next();
			String compartmentAlias = elm.getAttributeValue(PantherConstants.COMPARTMENT_ATT_NAME);
			GKInstance cmptInstance = compartAliasMap.get(compartmentAlias);
			String speciesID = elm.getAttributeValue(PantherConstants.SPECIES_ATT_NAME);
			GKInstance speciesInstance = idToInstanceMap.get(speciesID);
			if (speciesInstance == null)
				throw new IllegalStateException("PantherToReactomeConverter.handleActiveStates(): "
				      + "cannot find an instance for species id \"" + speciesID + "\".");
			String alias = elm.getAttributeValue(PantherConstants.ID_ATT_NAME);
			String activityText = elm.getChildTextTrim(PantherConstants.ACTIVITY_ELM_NAME,
			      PantherConstants.CELL_DESIGNER_NS);
			if (activityText == null || activityText.equals(PantherConstants.INACTIVE_STATE)) {
				// Default is inactivity
				idToInstanceMap.put(alias, speciesInstance);
			} else if (activityText.equals(PantherConstants.ACTIVE_STATE)) {
				// The same active state might be listed multiple times in a
				// pathway if a protein
				// is used in multiple reactions. See Akt in Angiogenesis.xml.
				GKInstance activeEntity = activeMap.get(speciesInstance);
				if (activeEntity == null) {
					// Need to clone speciesInstance
					activeEntity = instanceGenerator.cloneEntity(speciesInstance);
					activeMap.put(speciesInstance, activeEntity);
					List names = activeEntity.getAttributeValuesList(ReactomeJavaConstants.name);
					if (names != null && names.size() > 0) {
						String firstName = (String) names.get(0);
						firstName = "Active " + firstName;
						names.set(0, firstName); // Modify the first name
						InstanceDisplayNameGenerator.setDisplayName(activeEntity);
					}
				}
				idToInstanceMap.put(alias, activeEntity);
			}
		}
	}

	private Map<String, GKInstance> extractCompartmentAliases(Element modelElm) throws JaxenException {
		Map<String, GKInstance> rtn = new HashMap<String, GKInstance>();
		String xPath = ".//celldesigner:compartmentAlias";
		XPath myPath = XPathHelper.getXPathObject(xPath, modelElm);
		List elements = myPath.selectNodes(modelElm);
		if (elements == null || elements.size() == 0)
			return rtn;
		for (Iterator it = elements.iterator(); it.hasNext();) {
			Element elm = (Element) it.next();
			String alias = elm.getAttributeValue(PantherConstants.ID_ATT_NAME);
			String cmptID = elm.getAttributeValue(PantherConstants.COMPARTMENT_ATT_NAME);
			GKInstance cmptInstance = idToInstanceMap.get(cmptID);
			// No compartment Instance for default
			if (cmptInstance == null)
				continue;
			rtn.put(alias, cmptInstance);
		}
		return rtn;
	}

	private void processReactionsElement(Element reactionsElm) throws Exception {
		List reactionElmList = reactionsElm.getChildren(PantherConstants.REACTION_ELM_NAME,
		      PantherConstants.SBML_NS);
		if (reactionElmList == null || reactionElmList.size() == 0)
			return;
		Element rxtElm = null;
		List<Element> inhibitionElms = new ArrayList<Element>();
		GKInstance rxt = null;
		Map<String, GKInstance> reactionMap = new HashMap<String, GKInstance>();
		for (Iterator it = reactionElmList.iterator(); it.hasNext();) {
			rxtElm = (Element) it.next();
			if (isInhibition(rxtElm))
				inhibitionElms.add(rxtElm);
			else {
				String name = rxtElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
				if (name == null) {
					name = rxtElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
				}
				idToReactionMap.put(name, createReaction(rxtElm, reactionMap));
				// createReaction(rxtElm, reactionMap);
			}
		}
		for (Element elm : inhibitionElms) {
			handleSpeciesInhibition(elm);
		}
	}

	private void handleSpeciesInhibition(Element rxtElm) throws Exception {
		String id = rxtElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
		if (id == null)
			id = rxtElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
		List<GKInstance> inputInstances = extractInputs(rxtElm, id);
		List<GKInstance> outputInstances = extractOutputs(rxtElm, id);
		// Assume such an inhibition can occur in one-to-one mode
		if (inputInstances.size() > 1 || outputInstances.size() > 1)
			throw new IllegalStateException("PantherToReactomeConverter.handleSpeciesInhibition(): "
			      + "an inhibition has more than one input or output (" + id + ").");
		GKInstance inhibitor = inputInstances.get(0);
		GKInstance output = outputInstances.get(0);
		// This information might be lost if output is involved in a reaction has
		// more than one input. CellDesigner species information is saved in definition
		String iDef = (String) inhibitor.getAttributeValue(ReactomeJavaConstants.definition);
		String id1 = recoverID(iDef);
		String oDef = (String) output.getAttributeValue(ReactomeJavaConstants.definition);
		String id2 = recoverID(oDef);
		String name = idToInstanceMap.get(id1).getDisplayName() + " -| " + idToInstanceMap.get(id2).getDisplayName();
		List<GKInstance> inhibitedEvents = searchInhibitedReactions(output);
		if (inhibitedEvents.size() > 0) {
			for (GKInstance event : inhibitedEvents) {
				GKInstance nr = createInstance(ReactomeJavaConstants.NegativeRegulation);
				// RegulationType has been deleted.
//				GKInstance regulationType = instanceGenerator.getRegulationType(PantherConstants.INHIBITION_TYPE);
//				nr.setAttributeValue(ReactomeJavaConstants.regulationType, regulationType);
				nr.setAttributeValue(ReactomeJavaConstants.regulator, inhibitor);
				nr.setAttributeValue(ReactomeJavaConstants.regulatedEntity, event);
				nr.addAttributeValue(ReactomeJavaConstants.name, name);
				InstanceDisplayNameGenerator.setDisplayName(nr);
			}
		} else {
			// Sometimes the inhibited species doesn't participate any
			// reactions. See Akt
			// in Angiogensis.xml for an example. Wrap both input and output in
			// a reaction.
			GKInstance reaction = createInstance(ReactomeJavaConstants.Reaction);
			reaction.addAttributeValue(ReactomeJavaConstants.input, inhibitor);
			reaction.addAttributeValue(ReactomeJavaConstants.input, output);
			reaction.addAttributeValue(ReactomeJavaConstants.name, name);
			InstanceDisplayNameGenerator.setDisplayName(reaction);
			addReaction(reaction);
			String id2rname = rxtElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
			if (id2rname == null) {
				id2rname = rxtElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
			}
			idToReactionMap.put(id2rname, reaction);
		}
	}

	private String recoverID(String definition) {
		String id = null;
		int index = definition.indexOf(",");
		if (index < 0)
			id = definition;
		else
			id = definition.substring(0, index);
		return id;
	}

	private boolean isInhibition(Element rxtElm) throws JaxenException {
		String xPath = ".//celldesigner:reactionType";
		XPath myPath = XPathHelper.getXPathObject(xPath, rxtElm);
		Element elm = (Element) myPath.selectSingleNode(rxtElm);
		if (elm == null)
			return false;
		String reactionType = elm.getTextTrim();
		if (reactionType.equals(PantherConstants.INHIBITION_TYPE))
			return true;
		return false;
	}

	private String grepDbIds(List<GKInstance> instances) {
		StringBuilder builder = new StringBuilder();
		List<Long> ids = new ArrayList<Long>();
		for (GKInstance instance : instances) {
			Long id = instance.getDBID();
			ids.add(id);
		}
		Collections.sort(ids);
		for (Iterator<Long> it = ids.iterator(); it.hasNext();) {
			builder.append(it.next());
			if (it.hasNext())
				builder.append(",");
		}
		return builder.toString();
	}

	/**
	 * Create a reaction GKInstance based on the information in rxtElm. However, a reaction whose type is
	 * INHIBITION will NOT be created. They will be handled after all other types reactions are created.
	 * @param rxtElm
	 * @return
	 * @throws Exception
	 */
	private GKInstance createReaction(Element rxtElm, Map<String, GKInstance> reactionMap) throws Exception {
		// id is used to construct a definition
		String id = rxtElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
		if (id == null)
			id = rxtElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
		// Names used in Panther are meaningless.
		List<GKInstance> inputInstances = extractInputs(rxtElm, id);
		List<GKInstance> outputInstances = extractOutputs(rxtElm, id);
		// Create a key
		String inputIds = grepDbIds(inputInstances);
		String outputIds = grepDbIds(outputInstances);
		String key = inputIds + "->" + outputIds;
		GKInstance rxtInstance = reactionMap.get(key);
		if (rxtInstance == null) {
			rxtInstance = createInstance(ReactomeJavaConstants.Reaction);
			rxtInstance.setAttributeValue(ReactomeJavaConstants.input, inputInstances);
			rxtInstance.setAttributeValue(ReactomeJavaConstants.output, outputInstances);
			addReaction(rxtInstance);
			// displayname will be used by inhibition so it should be set here.
			String name = rxtElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
			if (name != null && !name.startsWith("re") && !(name.startsWith("r") && name.length() < 5)) {
				// Most likely generate by CellDesigner automatically
				rxtInstance.addAttributeValue(ReactomeJavaConstants.name, name);
			}
			InstanceDisplayNameGenerator.setDisplayName(rxtInstance);
			reactionMap.put(key, rxtInstance);
		}
		Element annotElm = rxtElm.getChild(PantherConstants.ANNOTATION_ELM_NAME, PantherConstants.SBML_NS);
		List list = annotElm.getChildren();
		List<GKInstance> modifierPathways = new ArrayList<GKInstance>();
		SchemaAttribute caEntityAtt = fileAdaptor.getSchema().getClassByName(
		      ReactomeJavaConstants.CatalystActivity).getAttribute(ReactomeJavaConstants.physicalEntity);
		String reactionType = null;
		boolean isModifierDone = false;
		for (Iterator it = list.iterator(); it.hasNext();) {
			Element child = (Element) it.next();
			String elmName = child.getName();
			if (elmName.equals(PantherConstants.REACTION_TYPE_ELM_NAME)) {
				reactionType = child.getTextTrim();
			} else if (elmName.equals(PantherConstants.LIST_OF_MODIFICATION_ELM_NAME)) {
				isModifierDone = true;
				List modificationElms = child.getChildren(PantherConstants.MODIFICATION_ELM_NAME,
				      PantherConstants.CELL_DESIGNER_NS);
				for (Iterator it1 = modificationElms.iterator(); it1.hasNext();) {
					Element modElm = (Element) it1.next();
					String alias = modElm.getAttributeValue(PantherConstants.ALIASES_ATT_NAME);
					String type = modElm.getAttributeValue(PantherConstants.TYPE_ATT_NAME);
					GKInstance speciesInstance = idToInstanceMap.get(alias);
					if (speciesInstance == null)
						throw new IllegalStateException("PantherToReactomeConverter.createReaction(): "
						      + "cannot find modifier \"" + alias + "\" for reaction \"" + id + "\".");
					// There are 12 types in total: state_transition, known_transition_omitted,
					// unknown_transition, catalysis, unknown_catalysis, inhibition, unknown_inhibition,
					// transport, transcriptional_activation, transcriptional_inhibition,
					// translational_actication, translational_inhibition. Catalysis and Unknown_catalysis
					// will be treated as CatalystAcvativity, inhibition, unknown_inhibition,
					// transcriptional_inhibition, translational_inhibition will be converted as
					// NegativeRegulation, and all others will be treated as PositiveReguation
					if (type.equals(PantherConstants.CATALYSIS_TYPE)
					      || type.equals(PantherConstants.UNKNOWN_CATALYSIS_TYPE)) {
						if (caEntityAtt.isValidValue(speciesInstance)) {
							GKInstance ca = getCatalystActivity(speciesInstance);
							if (rxtInstance != null)
								rxtInstance.addAttributeValue(ReactomeJavaConstants.catalystActivity, ca);
						}
					} else if (type.equals(PantherConstants.INHIBITION_TYPE)
					      || type.equals(PantherConstants.UNKNOWN_INHIBITION_TYPE)
					      || type.equals(PantherConstants.TRANSCRIPTIONAL_INHIBITION_TYPE)
					      || type.equals(PantherConstants.TRANSLATIONAL_INHIBITION_TYPE)) {
						createRegulation(speciesInstance, rxtInstance, ReactomeJavaConstants.NegativeRegulation);
					} else {
						createRegulation(speciesInstance, rxtInstance, ReactomeJavaConstants.PositiveRegulation);
					}
				}
			}
		}
		if (!isModifierDone) {
			Element listOfModifiersElm = rxtElm.getChild(PantherConstants.LIST_OF_MODIFIERS_ELM_NAME,
			      PantherConstants.SBML_NS);
			if (listOfModifiersElm != null) {
				List modifiersElm = listOfModifiersElm.getChildren(PantherConstants.modifierSpeciesReference,
				      PantherConstants.SBML_NS);
				if (modifiersElm != null) {
					for (Iterator it1 = modifiersElm.iterator(); it1.hasNext();) {
						Element modElm = (Element) it1.next();
						String species = modElm.getAttributeValue(PantherConstants.SPECIES_ATT_NAME);
						GKInstance speciesInstance = idToInstanceMap.get(species);
						if (speciesInstance == null)
							throw new IllegalStateException(
							      "PantherToReactomeConverter.createReaction(): cannot find modifier \"" + species
							            + "\" for reaction \"" + id + "\".");
						// No type defined: use Regulation for the time being.
						// However, regulation should not be used since it is abstract.
						createRegulation(speciesInstance, rxtInstance, ReactomeJavaConstants.Regulation);
					}
				}
			}
		}
		if (rxtInstance != null) {
			String def = (String) rxtInstance.getAttributeValue(ReactomeJavaConstants.definition);
			if (def == null)
				def = id + ", " + reactionType;
			else
				def = id + ", " + def;
			def = def + ", " + topPathway.getDisplayName();
			rxtInstance.setAttributeValue(ReactomeJavaConstants.definition, def);
			processEventNodeElement(rxtElm, rxtInstance);
		}
		return rxtInstance;
	}

	/**
	 * Prepares the celldesigner coordinates for every alias. It does not use species, because single species
	 * entries can be used more than once with different aliases. the result will be stored in the private
	 * class-variable coordinates. CellDesigner4 uses usualView in their coordinates. these relative
	 * coordinates are added to the species coordinates to acquire precise box coordinates
	 * @param annotationElm Annotation Element, child element of model
	 * @throws JDOMException XPath <code>selectNodes</code> has thrown an exception
	 */
	private void extractCoordinates(Element annotationElm) throws JaxenException {
		coordinates = new HashMap<String, Map<String, String>>();

		// getting the coordinates for the entire pathway via "modelDisplay"
		// coordinates of entire pathway will be saved like other coordinates, but under the
		// key of ReactomeConstants.PathwayDiagram
		String modelDisplayPath = ".//celldesigner:modelDisplay";
		XPath myPath = XPathHelper.getXPathObject(modelDisplayPath, annotationElm);
		List list = myPath.selectNodes(annotationElm);
		Map<String, String> coordSet = new HashMap<String, String>();
		if (list != null && list.size() > 0) {
			Element modelDisplayElm = (Element) list.get(0);
			String sizeX = modelDisplayElm.getAttributeValue(PantherConstants.SIZEX_ATT_NAME);
			String sizeY = modelDisplayElm.getAttributeValue(PantherConstants.SIZEY_ATT_NAME);

			coordSet.put(PantherConstants.W_ATT_NAME, sizeX);
			coordSet.put(PantherConstants.H_ATT_NAME, sizeY);
			coordSet.put(PantherConstants.SPECIES_ATT_NAME, ReactomeJavaConstants.PathwayDiagram);

			coordinates.put(ReactomeJavaConstants.PathwayDiagram, coordSet);
		}

		// continuing with rest of the coordinates
		String speciesAliasPath = ".//celldesigner:speciesAlias";
		myPath = XPathHelper.getXPathObject(speciesAliasPath, annotationElm);
		list = myPath.selectNodes(annotationElm);
		if (list != null && list.size() > 0) {
			String alias;
			String species;
			for (Iterator it = list.iterator(); it.hasNext();) {
				alias = null;
				Element speciesElm = (Element) it.next();
				alias = speciesElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
				species = speciesElm.getAttributeValue(PantherConstants.SPECIES_ATT_NAME);
				Element innerAliasesListElm = speciesElm.getChild(
				      PantherConstants.LIST_OF_INNER_ALIASES_ELM_NAME, PantherConstants.CELL_DESIGNER_NS);

				// handling complex coordinates
				if (innerAliasesListElm != null) {
					List<Element> innerAliases = innerAliasesListElm.getChildren(
					      PantherConstants.INNERALIAS_ELM_NAME, PantherConstants.CELL_DESIGNER_NS);
					for (Element aliasElm : innerAliases) {
						// these are the listOfInnerAliases Entries, innerAlias,
						// tagged inner0-innerx
						coordSet = new HashMap<String, String>();
						Box boundsBox = extractBox(aliasElm);
						String innerId = aliasElm.getAttributeValue(PantherConstants.HETERODIMERENTRY_ATT_NAME);
						coordSet.put(PantherConstants.SPECIES_ATT_NAME, getSpecies(species, innerId));
						coordSet.put(PantherConstants.X_ATT_NAME, boundsBox.getX().toString());
						coordSet.put(PantherConstants.Y_ATT_NAME, boundsBox.getY().toString());
						coordSet.put(PantherConstants.W_ATT_NAME, boundsBox.getW().toString());
						coordSet.put(PantherConstants.H_ATT_NAME, boundsBox.getH().toString());
						if ((alias != null || alias.length() == 0) && !coordSet.isEmpty()) {
							coordinates.put(alias + "," + innerId, coordSet);
						}
					}
				}
				// handling simple entity coordinates
				else {
					// these default species aliases, element name speciesAlias
					coordSet = new HashMap<String, String>();
					Box boundsBox = extractBox(speciesElm);
					Element boundsElm = speciesElm.getChild(PantherConstants.BOUNDS_ELM_NAME,
					      PantherConstants.CELL_DESIGNER_NS);
					String x = boundsElm.getAttributeValue(PantherConstants.X_ATT_NAME);
					String y = boundsElm.getAttributeValue(PantherConstants.Y_ATT_NAME);
					String width = boundsElm.getAttributeValue(PantherConstants.W_ATT_NAME);
					String height = boundsElm.getAttributeValue(PantherConstants.H_ATT_NAME);
					if (x != null && y != null && width != null && height != null) {
						x = doubleToIntString(x);
						y = doubleToIntString(y);
						width = doubleToIntString(width);
						height = doubleToIntString(height);

						coordSet.put(PantherConstants.SPECIES_ATT_NAME, species);
						coordSet.put(PantherConstants.X_ATT_NAME, x);
						coordSet.put(PantherConstants.Y_ATT_NAME, y);
						coordSet.put(PantherConstants.W_ATT_NAME, width);
						coordSet.put(PantherConstants.H_ATT_NAME, height);
						if ((alias != null || alias.length() == 0) && !coordSet.isEmpty()) {
							coordinates.put(alias, coordSet);
						}
					}
				}
			}
		}
		processInnerAliasCoordinates();
	}

	/**
	 * Takes inner aliases and summarises them to create a general alias with coordinates for a certain
	 * complex. inner aliases will remain, but generalised complex info will be added.
	 */
	private void processInnerAliasCoordinates() {
		Map<String, Map<String, String>> addCoordinates = new HashMap<String, Map<String, String>>();
		for (Iterator it = coordinates.keySet().iterator(); it.hasNext();) {
			String alias = (String) it.next();
			if (alias.contains("inner")) {
				String newAlias = alias.substring(0, alias.indexOf(","));
				if (addCoordinates.get(newAlias) == null) {
					String species = coordinates.get(alias).get("species");
					String newSpecies = species.substring(0, species.indexOf(","));
					Map<String, String> coordSet = new HashMap<String, String>();
					coordSet.put("species", newSpecies);
					addCoordinates.put(newAlias, coordSet);
				}
			}
		}
		coordinates.putAll(addCoordinates);
		for (Iterator it = coordinates.keySet().iterator(); it.hasNext();) {
			String alias = (String) it.next();
			if (!alias.contains("inner")) {
				continue;
			}
			String newAlias = alias.substring(0, alias.indexOf(","));
			Map<String, String> innerCoords = coordinates.get(alias);
			Map<String, String> coordSet = coordinates.get(newAlias);
			if (coordSet.size() <= 1) {
				coordSet.put(PantherConstants.X_ATT_NAME, innerCoords.get("x"));
				coordSet.put(PantherConstants.Y_ATT_NAME, innerCoords.get("y"));
				coordSet.put(PantherConstants.W_ATT_NAME, innerCoords.get("w"));
				coordSet.put(PantherConstants.H_ATT_NAME, innerCoords.get("h"));
			} else {
				for (Iterator it2 = innerCoords.keySet().iterator(); it2.hasNext();) {
					String key = (String) it2.next();
					if (!key.equals("species")) {
						Integer value = Integer.parseInt(innerCoords.get(key));
						if (key.equals("w") || key.equals("h")) {
							if (Integer.parseInt(coordSet.get(key)) < value) {
								coordSet.put(key, value.toString());
							}
						} else {
							if (Integer.parseInt(coordSet.get(key)) > value) {
								coordSet.put(key, value.toString());
							}
						}

					}
				}
			}
			coordinates.put(newAlias, coordSet);
		}
	}

	/**
	 * Small, rather futile string concat subfunction
	 */
	private String getSpecies(String species, String innerId) {
		return species + "," + innerId;
	}

	/**
	 * Extracts x, y, width and height from a celldesigner element and returns the data wrapped in an
	 * org.reactome.panther.Box-Class.
	 * @param elm Element to extract data from
	 * @return Box that contains x, y, width and height
	 */
	private Box extractBox(Element elm) {
		Element boundsElm = elm.getChild(PantherConstants.BOUNDS_ELM_NAME, PantherConstants.CELL_DESIGNER_NS);
		if (boundsElm != null) {
			Box boundsBox = new Box();
			boundsBox.setX(doubleStringToInt(boundsElm.getAttributeValue(PantherConstants.X_ATT_NAME)));
			boundsBox.setY(doubleStringToInt(boundsElm.getAttributeValue(PantherConstants.Y_ATT_NAME)));
			boundsBox.setW(doubleStringToInt(boundsElm.getAttributeValue(PantherConstants.W_ATT_NAME)));
			boundsBox.setH(doubleStringToInt(boundsElm.getAttributeValue(PantherConstants.H_ATT_NAME)));
			return boundsBox;
		}
		return null;
	}

	/**
	 * transforms and rounds a string in double format to a string in integer format
	 * @param doubleValue value to be transformed
	 * @return transformed integer value
	 */
	private String doubleToIntString(String doubleValue) {
		Double value = Double.parseDouble(doubleValue);
		Integer intValue = Math.round(value.floatValue());
		return intValue.toString();
	}

	/**
	 * transforms and rounds a string in double format to an integer
	 * @param doubleValue value to be transformed
	 * @return transformed integer value
	 */
	private int doubleStringToInt(String doubleValue) {
		Double value = Double.parseDouble(doubleValue);
		Integer intValue = Math.round(value.floatValue());
		return intValue;
	}

	/**
	 * Creates a PathwayDiagram GKInstance object.
	 * @param pathway - the GKInstance of the according pathway
	 * @return the PathwayDiagram GKInstance
	 * @throws Exception GKInstance-related Exception
	 */
	private GKInstance createPathwayDiagram(GKInstance pathway) throws Exception {
		String clsType = ReactomeJavaConstants.PathwayDiagram;
		GKInstance pathwayDiagram = createInstance(clsType);
		if (pathwayDiagram.getSchemClass().isValidAttribute(ReactomeJavaConstants.representedPathway)) {
			pathwayDiagram.addAttributeValue(ReactomeJavaConstants.representedPathway, pathway);
		}
		pathwayDiagram.addAttributeValue(ReactomeJavaConstants._displayName, "pathwayDiagram for "
		      + pathway.getAttributeValue(ReactomeJavaConstants._displayName));
		String width = coordinates.get(clsType).get(PantherConstants.W_ATT_NAME);
		String height = coordinates.get(clsType).get(PantherConstants.H_ATT_NAME);
		pathwayDiagram.addAttributeValue(ReactomeJavaConstants.width, Integer.parseInt(width));
		pathwayDiagram.addAttributeValue(ReactomeJavaConstants.height, Integer.parseInt(height));
		return pathwayDiagram;
	}

	/**
	 * extracts aliases and puts them into <code>reactionToAliasMap</code>
	 * @param reactionsElm listOfReactions DOM Element
	 * @throws JDOMException Error when parsing the document
	 */
	private void extractAliasMatching(Element reactionsElm) throws JaxenException {
		List list = reactionsElm.getChildren(PantherConstants.REACTION_ELM_NAME, PantherConstants.SBML_NS);
		if (list != null && list.size() > 0) {
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element reactionElm = (Element) it.next();
				String reactionID = reactionElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
				if (reactionID == null) {
					reactionID = reactionElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
				}
				Element annElm = reactionElm.getChild(PantherConstants.ANNOTATION_ELM_NAME,
				      PantherConstants.SBML_NS);

				List<String> baseInput = getBaseAliases(reactionElm, PantherConstants.BASE_REACTANTS_ELM_NAME,
				      PantherConstants.LIST_OF_REACTANTS_ELM_NAME);
				List<String> baseOutput = getBaseAliases(reactionElm, PantherConstants.BASE_PRODUCTS_ELM_NAME,
				      PantherConstants.LIST_OF_PRODUCTS_ELM_NAME);

				String inputAliasPath = ".//celldesigner:listOfReactantLinks/celldesigner:reactantLink";
				String outputAliasPath = ".//celldesigner:listOfProductLinks/celldesigner:productLink";
				String catalystPath = ".//celldesigner:listOfModification/celldesigner:modification";

				List<String> input = getAliases(annElm, inputAliasPath);
				List<String> output = getAliases(annElm, outputAliasPath);
				List<String> catalyst = getAliases(annElm, catalystPath);

				Map<String, List<String>> reaction = new HashMap<String, List<String>>();
				reaction.put("baseInput", baseInput);
				reaction.put("baseOutput", baseOutput);
				reaction.put("input", input);
				reaction.put("output", output);
				reaction.put("catalyst", catalyst);

				reactionToAliasMap.put(reactionID, reaction);
			}
		}
	}

	/**
	 * Extracts aliases of an element according to a certain xPath
	 * @param annotationElm
	 * @param xPath
	 * @return
	 */
	private List<String> getAliases(Element annotationElm, String xPath) {
		XPath myPath = null;
		List list = null;
		try {
			myPath = XPathHelper.getXPathObject(xPath, annotationElm);
			list = myPath.selectNodes(annotationElm);
		} catch (JaxenException ex) {
			ex.printStackTrace();
		}
		List<String> out = new ArrayList<String>();
		if (list != null && list.size() != 0) {
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element elm = (Element) it.next();
				String alias = getAliasAttributeValue(elm);
				out.add(alias);
			}
		}

		return out;
	}

	/**
	 * Checks an element and returns if there is an attribute called ALIAS or ALIASES respectively
	 * @param elm Element - The Element in Question
	 * @return The Attribute Value
	 */
	private String getAliasAttributeValue(Element elm) {
		List attributeList = elm.getAttributes();
		for (Iterator it = attributeList.iterator(); it.hasNext();) {
			Attribute attribute = (Attribute) it.next();
			if (attribute.getName().equals(PantherConstants.ALIASES_ATT_NAME)
			      || attribute.getName().equals(PantherConstants.ALIAS_ATT_NAME)) {
				return attribute.getValue();
			}
		}
		return "";
	}

	/**
	 * Retreives the base Alias of a reactant or a product
	 * @param reactionElm Reaction-Element
	 * @param listName Panther Specific element name for the reactantlist or productlist
	 * @return List of Aliases
	 */
	private List<String> getBaseAliases(Element reactionElm, String baseName, String listName) {
		Element annElm = reactionElm.getChild(PantherConstants.ANNOTATION_ELM_NAME, PantherConstants.SBML_NS);
		String inputString = annElm.getChildText(baseName, PantherConstants.CELL_DESIGNER_NS);
		List<String> speciesList = separate(inputString, ",");
		List<String> aliases = new ArrayList<String>();
		for (String species : speciesList) {
			String xPath = ".//" + listName + "/speciesReference[@species='" + species
			      + "']/annotation/celldesigner:alias";
			XPath myPath = null;
			List list = null;
			try {
				myPath = XPathHelper.getXPathObject(xPath, reactionElm);
				list = myPath.selectNodes(reactionElm);
			} catch (JaxenException e) {
				e.printStackTrace();
			}
			if (list != null && list.size() > 0) {
				Element annotation = (Element) list.get(0);
				aliases.add(annotation.getText());
			}
		}
		return aliases;
	}

	/**
	 * Separates a given String <code>toSeparate</code> where it matches <code>separator</code> and returns
	 * the result in list-form
	 * @param toSeparate - String
	 * @param separator - String
	 * @return List&lt;String> - List of separated Strings, without separator
	 */
	private List<String> separate(String toSeparate, String separator) {
		String rest = toSeparate;
		List<String> out = new ArrayList<String>();
		while (rest.indexOf(separator) != -1) {
			out.add(rest.substring(0, rest.indexOf(separator)));
			rest = rest.substring(rest.indexOf(separator) + 1);
		}
		out.add(rest);
		return out;
	}

	/**
	 * Generates vertices and edges for the current pathway
	 * @param pathwayDiagram the GKInstance for the pathway Diagram used here
	 */
	private void generateVertices(GKInstance pathwayDiagram) throws Exception {
		Set<String> reactionIds = idToReactionMap.keySet();
		GKInstance reactionVertex = null;
		for (String tmp : reactionIds) {
			if (tmp!=null) {
				GKInstance reaction = idToReactionMap.get(tmp);
				Map<String, List<String>> reactionMap = reactionToAliasMap.get(tmp);
	
				reactionVertex = createInstance(ReactomeJavaConstants.ReactionVertex);
				reactionVertex.addAttributeValue(ReactomeJavaConstants.representedInstance, reaction);
				reactionVertex.addAttributeValue(ReactomeJavaConstants.pathwayDiagram, pathwayDiagram);
	
				Line reactionCoords = handleVerticesAndChains(reactionMap, reactionVertex, pathwayDiagram);
				reactionVertex.addAttributeValue(ReactomeJavaConstants.x, reactionCoords.getAvgX());
				reactionVertex.addAttributeValue(ReactomeJavaConstants.y, reactionCoords.getAvgY());
	
				// reactions in Imre's Model all have width/height 8, Does that have to be changed?
				// Integer width = Integer.parseInt(coordinates.get(tmp).get(PantherConstants.W_ATT_NAME));
				// Integer height = Integer.parseInt(coordinates.get(tmp).get(PantherConstants.H_ATT_NAME));
				reactionVertex.addAttributeValue(ReactomeJavaConstants.width, 8);
				reactionVertex.addAttributeValue(ReactomeJavaConstants.height, 8);
			}
		}
	}

	/**
	 * Creates a vertex and connects it to <code>reactionVertex</code> in a new <code>Edge</code>, then
	 * returns the coordinates it used in a <code>Line</code>
	 * @param reactionMap A map containing the reactions with base and normal inputs and outputs, as well as
	 *           the catalyzer
	 * @param reactionVertex GKInstance for the reactionVertex that was already created
	 * @param pathwayDiagram GKInstance for the pathway Diagram of the reaction
	 * @return The calculated line coordinates for the reaction in question
	 * @throws Exception
	 */
	private Line handleVerticesAndChains(Map<String, List<String>> reactionMap, GKInstance reactionVertex,
	      GKInstance pathwayDiagram) throws Exception {
		Set<String> keys = reactionMap.keySet();
		Line reactionCoords = new Line();
		for (String type : keys) {
			List<String> typeList = reactionMap.get(type);
			for (Iterator it = typeList.iterator(); it.hasNext();) {
				String alias = (String) it.next();
				GKInstance entity = idToInstanceMap.get(alias);
				Map<String, String> coordSet = coordinates.get(alias);
				Integer x = Integer.parseInt(coordSet.get(PantherConstants.X_ATT_NAME));
				Integer y = Integer.parseInt(coordSet.get(PantherConstants.Y_ATT_NAME));
				Integer w = Integer.parseInt(coordSet.get(PantherConstants.W_ATT_NAME));
				Integer h = Integer.parseInt(coordSet.get(PantherConstants.H_ATT_NAME));

				if (type.contains("baseInput")) {
					reactionCoords.addToSource(x, y);
				}
				if (type.equals("baseOutput")) {
					reactionCoords.addToTarget(x, y);
				}

				GKInstance entityVertex = createVertex(ReactomeJavaConstants.EntityVertex, entity,
				      pathwayDiagram, x, y, w, h);

				if (!type.equals("output") && !type.equals("baseOutput"))
					createEdgeItem(entityVertex, reactionVertex);
				else
					createEdgeItem(reactionVertex, entityVertex);
			}
		}
		return reactionCoords;
	}

	/**
	 * Creates a Vertex element with the given parameters
	 * @param classType decides if it will be a ReactionVertex or an EntityVertex
	 * @param representedInstance the GKInstance this vertex is representing
	 * @param pathwayDiagram GKInstance of the pathway Diagram this Vertex is in
	 * @param x x-Coordinate of the Vertex
	 * @param y y-Coordinate of the Vertex
	 * @param width width of the Vertex
	 * @param height height of the vertex
	 * @return the Vertex that was created
	 * @throws Exception Instance Exception
	 */
	private GKInstance createVertex(String classType, GKInstance representedInstance,
	      GKInstance pathwayDiagram, Integer x, Integer y, Integer width, Integer height) throws Exception {
		if (!instanceGenerator.isVertexType(classType))
			return null;
		GKInstance vertex = createInstance(classType);
		vertex.addAttributeValue(ReactomeJavaConstants.representedInstance, representedInstance);
		vertex.addAttributeValue(ReactomeJavaConstants.pathwayDiagram, pathwayDiagram);
		vertex.addAttributeValue(ReactomeJavaConstants.x, x);
		vertex.addAttributeValue(ReactomeJavaConstants.y, y);
		vertex.addAttributeValue(ReactomeJavaConstants.width, width);
		vertex.addAttributeValue(ReactomeJavaConstants.height, height);

		return vertex;
	}

	/**
	 * Creates a GKInstance of the type <code>Edge</code> from <code>sourceVertex</code> to
	 * <code>targetVertex</code>
	 * @param sourceVertex
	 * @param targetVertex
	 * @return Created Edge Instance or <code>null</code> if parameters are not of the type
	 *         <code>Vertex</code>
	 * @throws Exception Instance Exception
	 */
	private GKInstance createEdgeItem(GKInstance sourceVertex, GKInstance targetVertex) throws Exception {
		GKInstance edgeItem = null;
		if (sourceVertex.getSchemClass().isa("Vertex") && targetVertex.getSchemClass().isa("Vertex")) {
			edgeItem = createInstance(ReactomeJavaConstants.Edge);
			edgeItem.addAttributeValue(ReactomeJavaConstants.sourceVertex, sourceVertex);
			edgeItem.addAttributeValue(ReactomeJavaConstants.targetVertex, targetVertex);
		}
		return edgeItem;
	}

	/**
	 * @param rxtElm
	 * @param id
	 * @return
	 * @throws JaxenException
	 */
	private List<GKInstance> extractInputs(Element rxtElm, String id) throws JaxenException {
		// Handle reactants
		Element reactantsElm = rxtElm.getChild(PantherConstants.LIST_OF_REACTANTS_ELM_NAME,
		      PantherConstants.SBML_NS);
		String xPath = ".//celldesigner:alias";
		XPath myPath = XPathHelper.getXPathObject(xPath, reactantsElm);
		List list = myPath.selectNodes(reactantsElm);
		if (list != null && list.size() > 0) {
			List<GKInstance> rtn = new ArrayList<GKInstance>(list.size());
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element elm = (Element) it.next();
				String alias = elm.getTextTrim();
				GKInstance input = idToInstanceMap.get(alias);
				if (input == null)
					throw new IllegalStateException("PantherToReactomeConverter.createReaction(): "
					      + "cannot find input \"" + alias + "\" for reaction \"" + id + "\".");
				rtn.add(input);
			}
			return rtn;
		} else {
			// Try the convertional way
			list = reactantsElm.getChildren(PantherConstants.SPECIES_REFERENCE, PantherConstants.SBML_NS);
			List<GKInstance> rtn = new ArrayList<GKInstance>(list.size());
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element elm = (Element) it.next();
				String name = elm.getAttributeValue(PantherConstants.SPECIES_ATT_NAME);
				GKInstance input = idToInstanceMap.get(name);
				if (input == null)
					throw new IllegalStateException("PantherToReactomeConverter.createReaction(): "
					      + "cannot find input \"" + name + "\" for reaction \"" + id + "\".");
				rtn.add(input);
			}
			return rtn;
		}
	}

	/**
	 * TODO: this method should be merged with extractInputs().
	 * 
	 * @param rxtElm
	 * @param id
	 * @return
	 * @throws JaxenException
	 */
	private List<GKInstance> extractOutputs(Element rxtElm, String id) throws JaxenException {
		// Handle reactants
		Element productsElm = rxtElm.getChild(PantherConstants.LIST_OF_PRODUCTS_ELM_NAME,
		      PantherConstants.SBML_NS);
		String xPath = ".//celldesigner:alias";
		XPath myPath = XPathHelper.getXPathObject(xPath, productsElm);
		List list = myPath.selectNodes(productsElm);
		if (list != null && list.size() > 0) {
			List<GKInstance> rtn = new ArrayList<GKInstance>(list.size());
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element elm = (Element) it.next();
				String alias = elm.getTextTrim();
				GKInstance output = idToInstanceMap.get(alias);
				if (output == null)
					throw new IllegalStateException("PantherToReactomeConverter.createReaction(): "
					      + "cannot find output \"" + alias + "\" for reaction \"" + id + "\".");
				rtn.add(output);
			}
			return rtn;
		} else {
			// Try the convertional way
			list = productsElm.getChildren(PantherConstants.SPECIES_REFERENCE, PantherConstants.SBML_NS);
			List<GKInstance> rtn = new ArrayList<GKInstance>(list.size());
			for (Iterator it = list.iterator(); it.hasNext();) {
				Element elm = (Element) it.next();
				String name = elm.getAttributeValue(PantherConstants.SPECIES_ATT_NAME);
				GKInstance output = idToInstanceMap.get(name);
				if (output == null)
					throw new IllegalStateException("PantherToReactomeConverter.createReaction(): "
					      + "cannot find output \"" + name + "\" for reaction \"" + id + "\".");
				rtn.add(output);
			}
			return rtn;
		}
	}

	private void createRegulation(GKInstance regulator, GKInstance target, String clsType) throws Exception {
		if (target == null)
			return; // Just in case
		GKInstance regulation = createInstance(clsType);
		regulation.addAttributeValue(ReactomeJavaConstants.regulatedEntity, target);
		regulation.addAttributeValue(ReactomeJavaConstants.regulator, regulator);
		InstanceDisplayNameGenerator.setDisplayName(regulation);
	}

	private GKInstance getCatalystActivity(GKInstance catalyst) throws Exception {
		GKInstance ca = caMap.get(catalyst);
		if (ca == null) {
			ca = createInstance(ReactomeJavaConstants.CatalystActivity);
			ca.addAttributeValue(ReactomeJavaConstants.physicalEntity, catalyst);
			InstanceDisplayNameGenerator.setDisplayName(ca);
			caMap.put(catalyst, ca);
		}
		return ca;
	}

	private Map<String, Protein> extractProteins(Element modelElm) throws JaxenException {
		Map<String, Protein> proteinMap = new HashMap<String, Protein>();
		String xPath = ".//celldesigner:protein";
		Element annotationElm = modelElm.getChild(PantherConstants.ANNOTATION_ELM_NAME,
		      PantherConstants.SBML_NS);
		XPath myPath = XPathHelper.getXPathObject(xPath, annotationElm);
		List proteinElms = myPath.selectNodes(annotationElm);
		Element elm = null;
		for (Iterator it = proteinElms.iterator(); it.hasNext();) {
			elm = (Element) it.next();
			String proteinRef = elm.getAttributeValue(PantherConstants.ID_ATT_NAME);
			String name = elm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
			name = nameHandler.cleanupName(name);
			String type = elm.getAttributeValue(PantherConstants.TYPE_ATT_NAME);
			Protein protein = new Protein(proteinRef, name, type);
			proteinMap.put(proteinRef, protein);
			// Check if there are any positions defined
			Element modificationList = elm.getChild(PantherConstants.LIST_OF_MODIFICATION_RESIDUE_ELM_NAME,
			      PantherConstants.CELL_DESIGNER_NS);
			if (modificationList == null)
				continue;
			List list = modificationList.getChildren(PantherConstants.MODIFICATION_RESIDUE_ELM_NAME,
			      PantherConstants.CELL_DESIGNER_NS);
			if (list == null && list.size() == 0)
				continue;
			for (Iterator it1 = list.iterator(); it1.hasNext();) {
				Element tmp = (Element) it1.next();
				String id = tmp.getAttributeValue(PantherConstants.ID_ATT_NAME);
				name = tmp.getAttributeValue(PantherConstants.NAME_ATT_NAME);
				// name might be null
				protein.addPositions(id, name);
			}
		}
		return proteinMap;
	}

	private void processSpeciesElement(Element speciesListElm) throws Exception {
		List elmList = speciesListElm.getChildren(PantherConstants.SPECIES_ELM_NAME, PantherConstants.SBML_NS);
		if (elmList == null || elmList.size() == 0)
			return; // Should NOT be a case
		Element elm = null;
		List<Element> complexElms = new ArrayList<Element>();
		// The first loop run to handle non-complex species
		for (Iterator it = elmList.iterator(); it.hasNext();) {
			elm = (Element) it.next();
			GKInstance instance = createEntity(elm);
			// idToEntityMap.put(elm.getAttributeValue(PantherConstants.ID_ATT_NAME),
			// instance);
			if (instance.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
				complexElms.add(elm); // It must be an complex
			}
		}
		for (Element complexElm : complexElms) {
			processComplexComponents(complexElm);
		}
		List<GKInstance> complexList = (List<GKInstance>) fileAdaptor
		      .fetchInstancesByClass(ReactomeJavaConstants.Complex);
		for (GKInstance instance : complexList) {
			String name = (String) instance.getAttributeValue("name");
			String compartment = "";
			String complexName = "";
			if (name.startsWith("s") && name.length() < 5) {
				try {
					Integer.parseInt(name.substring(1));
					compartment = instance.getDisplayName().substring(name.length());
					List<GKInstance> members = instance.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
					for (GKInstance memberInstance : members) {
						if (complexName.length() != 0)
							complexName += ":";
						complexName += memberInstance.getAttributeValue("name").toString();
					}
				} catch (NumberFormatException ex) {
				}
			} else {
				complexName = nameHandler.cleanupComplexName(name);
			}
			instance.setAttributeValue(ReactomeJavaConstants.name, complexName);
			instance.setDisplayName(complexName + compartment);
		}
	}

	private void processComplexComponents(Element complexElm) throws Exception {
		String id = complexElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
		if (id == null)
			id = complexElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
		GKInstance complex = idToInstanceMap.get(id);
		String xPath = ".//celldesigner:" + PantherConstants.HETERODIMER_ENTRY_ELM_NAME;
		XPath myPath = XPathHelper.getXPathObject(xPath, complexElm);
		List entries = myPath.selectNodes(complexElm);
		if (entries == null || entries.size() == 0)
			return;
		for (Iterator it = entries.iterator(); it.hasNext();) {
			Element entryElm = (Element) it.next();
			handleComplexComp(entryElm, complex);
		}
		// If a complex component is not defined as a species, its Panther ID
		// should
		// be extracted from Heterodimer Member Info in note element
		Element noteElm = complexElm.getChild(PantherConstants.NOTES_ELM_NAME, complexElm.getNamespace());
		String[] text = getTextFromNoteElm(noteElm);
		// Check for Heterodimer Member Info
		processHeterodimerMemberInfo(text, complex);
	}

	private void processHeterodimerMemberInfo(String[] text, GKInstance complex) throws Exception {
		if (text == null || text.length == 0) {
			System.out.println("Complex has no note: " + complex);
			return;
		}
		String header = "Heterodimer Member Info";
		List components = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
		for (String line : text) {
			if (line.startsWith(header)) {
				String sub = line.substring(header.length() + 2).trim(); // Escape ":"
				String[] tokens = sub.split(";"); // each token for each component
				boolean found = false;
				for (String token : tokens) {
					String[] info = token.split("#"); // Work with four elements only
					if (info.length < 4)
						continue;
					// First is name, last is the Panther id Hard code this to make sure it is correct
					String name = info[0];
					// In Heterotrimeric_G_protein_signaling_pathway_Gi_and_Gs_mediated_pathway.xml
					// there are no _ between alpha, beta, and gamma
					String id = info[3];
					found = false;
					// Find component for name
					for (Iterator it = components.iterator(); it.hasNext();) {
						GKInstance comp = (GKInstance) it.next();
						List names = comp.getAttributeValuesList(ReactomeJavaConstants.name);
						if (names.contains(name)) {
							found = true;
							GKInstance dbRef = (GKInstance) comp
							      .getAttributeValue(ReactomeJavaConstants.crossReference);
							if (dbRef == null) {
								// Create one
								dbRef = instanceGenerator.getDBIdentifier(id);
								comp.setAttributeValue(ReactomeJavaConstants.crossReference, dbRef);
							}
							break;
						}
					}
					if (!found) {
						// Check if crossReference is needed
						for (Iterator it = components.iterator(); it.hasNext();) {
							GKInstance comp = (GKInstance) it.next();
							if (comp.getAttributeValue(ReactomeJavaConstants.crossReference) == null) {
								System.out.println("Component " + name + " in complex " + complex.getDBID()
								      + " cannot be mapped!");
								break;
							}
						}
					}
				}
				break;
			}
		}
	}

	private void handleComplexComp(Element entryElm, GKInstance complex) throws Exception {
		List children = entryElm.getChildren();
		Element childElm = null;
		String clsName = null;
		String name = null;
		String proteinRef = null;
		String stoiStr = null;
		Element modificationElm = null;
		for (Iterator it = children.iterator(); it.hasNext();) {
			childElm = (Element) it.next();
			String elmName = childElm.getName();
			if (elmName.equals(PantherConstants.CLASS_ELM_NAME)) {
				clsName = childElm.getTextTrim();
			} else if (elmName.equals(PantherConstants.PROTEIN_REFERENCE_ELM_NAME)) {
				proteinRef = childElm.getTextTrim();
			} else if (elmName.equals(PantherConstants.NAME_ELM_NAME)) {
				name = childElm.getTextTrim();
				name = nameHandler.cleanupName(name);
			} else if (elmName.equals(PantherConstants.STATE_ELM_NAME)) {
				stoiStr = childElm.getChildTextTrim(PantherConstants.HOMODIMER_ELM_NAME,
				      PantherConstants.CELL_DESIGNER_NS);
				modificationElm = childElm.getChild("listOfModifications", PantherConstants.CELL_DESIGNER_NS);
			}
		}
		GKInstance comp = null;
		GKInstance compartment = (GKInstance) complex.getAttributeValue(ReactomeJavaConstants.compartment);
		if (clsName.equals(PantherConstants.PROTEIN_TYPE_NAME)) {
			// proteinReference is used for PROTEIN class
			Protein protein = proteinMap.get(proteinRef);
			// Modifications should be handled too. Only Proteins can have
			// modifications
			List<GKInstance> modifications = null;
			if (modificationElm != null)
				modifications = getModifications(modificationElm, protein);
			List<GKInstance> entities = protein.getEntities();
			if (entities == null || entities.size() == 0) {
				comp = instanceGenerator.createEntity(protein, topPathway.getDisplayName());
				if (compartment != null) {
					comp.addAttributeValue(ReactomeJavaConstants.compartment, compartment);
				}
				if (modifications != null)
					comp.setAttributeValue(ReactomeJavaConstants.hasModifiedResidue, modifications);
				InstanceDisplayNameGenerator.setDisplayName(comp);
			} else {
				for (GKInstance entity : entities) {
					if (checkEntity(entity, compartment, modifications)) {
						comp = entity;
						break;
					}
				}
				if (comp == null) {
					comp = instanceGenerator.createEntity(protein, topPathway.getDisplayName());
					// want to copy DatabaseCrossReference
					GKInstance speciesInstance = entities.get(0);
					GKInstance dbcr = (GKInstance) speciesInstance
					      .getAttributeValue(ReactomeJavaConstants.crossReference);
					if (dbcr != null)
						comp.addAttributeValue(ReactomeJavaConstants.crossReference, dbcr);
					if (compartment != null)
						comp.addAttributeValue(ReactomeJavaConstants.compartment, compartment);
					if (modifications != null)
						comp.setAttributeValue(ReactomeJavaConstants.hasModifiedResidue, modifications);
					List names = speciesInstance.getAttributeValuesList(ReactomeJavaConstants.name);
					comp.setAttributeValue(ReactomeJavaConstants.name, new ArrayList(names));
					InstanceDisplayNameGenerator.setDisplayName(comp);
				}
			}
		} else {
			comp = nameToInstanceMap.get(name);
			if (comp == null) {
				comp = instanceGenerator.createEntity(clsName);
				comp.addAttributeValue(ReactomeJavaConstants.name, name);
				if (compartment != null)
					comp.addAttributeValue(ReactomeJavaConstants.compartment, compartment);
				InstanceDisplayNameGenerator.setDisplayName(comp);
				nameToInstanceMap.put(name, comp);
			}
		}
		if (comp == null)
			throw new IllegalStateException(
			      "PantherToReactomeConverter.handleComplexComp(): cannot find complex component! ("
			            + (proteinRef != null ? proteinRef : name) + ")");
		SchemaAttribute hasCompAtt = complex.getSchemClass().getAttribute(ReactomeJavaConstants.hasComponent);
		if (!hasCompAtt.isValidValue(comp))
			return; // It might be a pathway
		if (stoiStr == null)
			complex.addAttributeValue(hasCompAtt, comp);
		else {
			int stoi = Integer.parseInt(stoiStr);
			for (int i = 0; i < stoi; i++)
				complex.addAttributeValue(hasCompAtt, comp);
		}
	}

	/**
	 * This helper method is used to check if a PhysicalEntity can be used as a component for a Complex. This
	 * is refactored from method handleComplexComp().
	 */
	private boolean checkEntity(GKInstance entity, GKInstance compartment, List<GKInstance> modifications)
	      throws Exception {
		GKInstance value = (GKInstance) entity.getAttributeValue(ReactomeJavaConstants.compartment);
		if (value != compartment)
			return false;
		List valueList = entity.getAttributeValuesList(ReactomeJavaConstants.hasModifiedResidue);
		// Make compare simple
		Set set1 = null;
		if (modifications == null)
			set1 = new HashSet();
		else
			set1 = new HashSet(modifications);
		Set set2 = null;
		if (valueList == null)
			set2 = new HashSet();
		else
			set2 = new HashSet(valueList);
		return set1.equals(set2);
	}

	private GKInstance createEntity(Element elm) throws Exception {
		String id = nameHandler.getSpeciesID(elm);
		String compartment = elm.getAttributeValue(PantherConstants.COMPARTMENT_ATT_NAME);
		Element annotationElm = elm.getChild(PantherConstants.ANNOTATION_ELM_NAME, PantherConstants.SBML_NS);
		String positionToCompartment = null;
		String classType = null;
		String proteinRef = null;
		List children = annotationElm.getChildren();
		Element stateElm = null;
		// Don't want to handle complex in this loop run
		for (Iterator it = children.iterator(); it.hasNext();) {
			Element child = (Element) it.next();
			String childName = child.getName();
			if (childName.equals(PantherConstants.POSITION_TO_COMPARTMENT_ELM_NAME)) {
				positionToCompartment = child.getTextTrim();
			} else if (childName.equals(PantherConstants.SPECIES_IDENTITY_ELM_NAME)) {
				classType = child.getChildTextTrim(PantherConstants.CLASS_ELM_NAME,
				      PantherConstants.CELL_DESIGNER_NS);
				if (classType.equals(PantherConstants.PROTEIN_TYPE_NAME)) {
					// Need the proteinreference
					proteinRef = child.getChildTextTrim(PantherConstants.PROTEIN_REFERENCE_ELM_NAME,
					      PantherConstants.CELL_DESIGNER_NS);
				}
				stateElm = child.getChild(PantherConstants.STATE_ELM_NAME, PantherConstants.CELL_DESIGNER_NS);
			} else if (childName.equals(PantherConstants.HETERODIMER_IDENTITY_ELM_NAME)) {
				classType = PantherConstants.COMPLEX_TYPE_NAME;
			}
		}
		if (classType == null) {
			classType = PantherConstants.UNKNOWN_TYPE_NAME;
		}
		// Handle names
		List<String> names = nameHandler.getEntityNames(elm, proteinMap);
		// Handle compartment
		GKInstance compInstance = instanceGenerator.getCompartment(compartment, positionToCompartment,
		      idToInstanceMap);
		GKInstance instance = null;
		boolean isOld = true;
		if (instanceGenerator.isSimpleType(classType))
			instance = instanceGenerator.searchSimpleEntity(classType, names, compInstance);
		if (instance == null) {
			instance = instanceGenerator.createEntity(classType);
			isOld = false;
		}
		if (instance == null)
			throw new IllegalStateException("PantherToReactomeConverter.createEntity(): "
			      + "cannot generate an instance for a species.");
		// Track degraded entities to assign a meaningful name afterwards
		if (classType.equals(PantherConstants.DEGRADED_TYPE_NAME))
			degradedEntities.add(instance);
		Protein protein = null;
		if (proteinRef != null) {
			protein = proteinMap.get(proteinRef);
			protein.addEntity(instance);
		}
		if (!isOld) {
			// Handle definition
			String definition = null;
			if (classType.equals(PantherConstants.COMPLEX_TYPE_NAME))
				definition = id;
			else
				definition = id + ", " + classType;
			if (protein != null) {
				definition = definition + ", " + protein.getType();
			}
			if (instance.getSchemClass().isa(ReactomeJavaConstants.GenomeEncodedEntity))
				definition += (", " + topPathway.getDisplayName()); // For debugging
			instance.setAttributeValue(ReactomeJavaConstants.definition, definition);
			if (names.size() > 0)
				instance.setAttributeValue(ReactomeJavaConstants.name, names);
			if (compInstance != null)
				instance.addAttributeValue(ReactomeJavaConstants.compartment, compInstance);
			// Handle noteElm after registering to map since the first name
			// should be used for nameToInstanceMap. Long Name might move celldesginer:name to
			// the second position. Need to get properties from notes
			Element noteElm = elm.getChild(PantherConstants.NOTES_ELM_NAME, PantherConstants.SBML_NS);
			parseEntityProperties(noteElm, instance);
			if (stateElm != null) {
				List<GKInstance> modifications = getModifications(stateElm, protein);
				if (modifications != null)
					instance.setAttributeValue(ReactomeJavaConstants.hasModifiedResidue, modifications);
			}
			cleanupEntityNames(instance);
			InstanceDisplayNameGenerator.setDisplayName(instance);
		}
		if (!classType.equals(PantherConstants.PROTEIN_TYPE_NAME) &&
		// Sometime a complex name is the same as its component: e.g. GPI linker
		      // in Plasminogen activating cascade.xml.
		      !classType.equals(PantherConstants.COMPLEX_TYPE_NAME)) {
			// In case these instances are used in complexes.
			// Key for nameToInstanceMap should be celldesigner name
			String firstName = (String) instance.getAttributeValue(ReactomeJavaConstants.name);
			if (firstName != null) // No names for complexes in level1. They will not be needed.
				nameToInstanceMap.put(firstName, instance);
		}
		idToInstanceMap.put(id, instance);
		return instance;
	}

	/**
	 * Sometimes there will be multiple names for single instances. This method will try to filter out
	 * Panther-style names like s123 because these are rather useless. If these Panther names remain, the name
	 * generator might use this name as a display name
	 * @param instance instance in question
	 * @throws Exception GKInstance-Exception
	 */
	private void cleanupEntityNames(GKInstance instance) throws Exception {
		List<String> names = instance.getAttributeValuesList(ReactomeJavaConstants.name);
		String name = names.size() > 0 ? names.get(0) : "";
		while (names.size() > 1 && name.length() < 5 && name.startsWith("s")) {
			try {
				Integer.parseInt(name.substring(1));
				names.remove(0);
				name = names.get(0);
			} catch (NumberFormatException ex) {
				break;
			}
		}
	}

	private List<GKInstance> getModifications(Element modificationElm, Protein protein) throws Exception {
		String xPath = ".//celldesigner:modification";
		XPath myPath = XPathHelper.getXPathObject(xPath, modificationElm);
		List modificationsList = myPath.selectNodes(modificationElm);
		if (modificationsList == null || modificationsList.size() == 0)
			return null;
		List<GKInstance> rtn = new ArrayList<GKInstance>();
		Element elm = null;
		for (Iterator it = modificationsList.iterator(); it.hasNext();) {
			elm = (Element) it.next();
			GKInstance modificationInstance = instanceGenerator.getModificationInstance(elm, protein);
			rtn.add(modificationInstance);
		}
		return rtn;
	}

	private void parseEntityProperties(Element noteElm, GKInstance entity) throws Exception {
		String[] lines = getTextFromNoteElm(noteElm);
		if (lines == null)
			return;
		nameHandler.getNamesFromNotes(lines, entity);
		for (String line : lines) {
			line = line.trim();
			if (line.length() == 0)
				continue;
			if (line.startsWith(PantherConstants.ACCESSION_LABEL)) {
				String accession = line.substring(PantherConstants.ACCESSION_LABEL.length() + 1).trim();
				GKInstance dbInstance = instanceGenerator.getDBIdentifier(accession);
				entity.addAttributeValue(ReactomeJavaConstants.crossReference, dbInstance);
			}
		}
	}

	private void processCompartmentsElement(Element compartmentsElm) throws Exception {
		List elmList = compartmentsElm.getChildren(PantherConstants.COMPARTMENT_ELM_NAME,
		      PantherConstants.SBML_NS);
		if (elmList == null || elmList.size() == 0)
			return;
		Element elm = null;
		String id = null;
		// Initialize the compartments in the first loop run
		// Keep outside relationships in this map
		Map<String, String> outsideMap = new HashMap<String, String>();
		for (Iterator it = elmList.iterator(); it.hasNext();) {
			elm = (Element) it.next();
			GKInstance instance = instanceGenerator.createCompartment(elm, nameHandler);
			if (instance == null)
				continue; // Default
			id = nameHandler.getSpeciesID(elm);
			idToInstanceMap.put(id, instance);
			String outside = elm.getAttributeValue(PantherConstants.OUTSIDE_ATT_NAME);
			if (outside != null) { // It is an optional attribute
				outsideMap.put(id, outside);
			}
		}
		// Figure out outside relationships here
		Set<String> ids = outsideMap.keySet();
		for (String tmp : ids) {
			String outside = outsideMap.get(tmp);
			if (outside.equals(PantherConstants.DEFAULT_COMPARTMENT_NAME))
				continue;
			GKInstance child = idToInstanceMap.get(tmp);
			GKInstance parent = idToInstanceMap.get(outside);
			List parents = child.getAttributeValuesList(ReactomeJavaConstants.componentOf);
			if (parents == null || parents.size() == 0)
				child.addAttributeValue(ReactomeJavaConstants.componentOf, parent);
			else if (!parents.contains(parent))
				child.addAttributeValue(ReactomeJavaConstants.componentOf, parent);
		}
	}

	private String[] getTextFromNoteElm(Element noteElm) {
		if (noteElm == null)
			return null; // Nothing to process
		String text = null;
		Element htmlElm = noteElm.getChild(PantherConstants.HTML_ELM_NAME, PantherConstants.HTML_NS);
		if (htmlElm == null) {
			text = noteElm.getTextTrim();
			// In some pathways (e.g. TGF-beta, text is attached to note directly)
		} else {
			Element bodyElm = htmlElm.getChild(PantherConstants.BODY_ELM_NAME, PantherConstants.HTML_NS);
			if (bodyElm == null)
				return null;
			text = bodyElm.getTextTrim();
		}
		if (text.equals(PantherConstants.EMPTY_PATHWAY_NOTE_MARK))
			return null;
		String[] lines = text.split("(\\n|\\r)");
		return lines;
	}

	private void extractPathwayProp(Element modelElm, GKInstance rPathway) throws Exception {
		String id = modelElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
		// for level1 model, name should be used
		if (id == null)
			id = modelElm.getAttributeValue(PantherConstants.NAME_ATT_NAME);
		id = nameHandler.cleanupName(id);
		// model id should be used as the name for GKInstance
		rPathway.setDisplayName(id);
		rPathway.addAttributeValue(ReactomeJavaConstants.name, id);
		processEventNodeElement(modelElm, rPathway);
	}

	/**
	 * This helper method is used to process notes element for both Pathway and Reaction.
	 * @param eventElm
	 * @param rEvent
	 * @throws Exception
	 */
	private void processEventNodeElement(Element eventElm, GKInstance rEvent) throws Exception {
		// Analyze notes for the pathway
		Element noteElm = eventElm.getChild(PantherConstants.NOTES_ELM_NAME, PantherConstants.SBML_NS);
		String[] lines = getTextFromNoteElm(noteElm);
		if (lines == null)
			return;
		StringBuilder summationText = new StringBuilder();
		// Get the text summation
		for (String line : lines) {
			if (line.trim().length() == 0)
				continue;
			line = line.trim();
			if (line.startsWith(PantherConstants.WEB_SITE_LABEL)
			      || line.startsWith(PantherConstants.FREE_TEXT_LABEL)
			      || line.startsWith(PantherConstants.MEDLINE_LABEL)
			      || line.startsWith(PantherConstants.PMID_LABEL)) {
				// Some pubmed ids start with Medline
				GKInstance literatureReferenece = instanceGenerator.createLiteratureReference(line);
				if (literatureReferenece != null)
					rEvent.addAttributeValue(ReactomeJavaConstants.literatureReference, literatureReferenece);
			} else {
				if (summationText.length() > 0)
					summationText.append("\n");
				summationText.append(line);
			}
		}
		// Construct a Summation GKInstance for pathway
		GKInstance summation = instanceGenerator.createSummation(summationText.toString());
		if (summation != null)
			rEvent.addAttributeValue(ReactomeJavaConstants.summation, summation);
	}

	private GKInstance createInstance(String clsName) throws Exception {
		return instanceGenerator.createInstance(clsName);
	}

	public void save(String projectFileName) throws Exception {
		this.getFileAdaptor().save(projectFileName);
	}
}
