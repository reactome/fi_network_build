/*
 * Created on Feb 28, 2006
 *
 */
package org.reactome.panther;

import org.jdom.Namespace;

/**
 * A setup class for elements and attributes in a Panther XML/SBML file 
 *
 */
public class PantherConstants {
    // A special name for Panther Database
    public static final String PANTHER_DB_NAME = "pantherdb";
    public static final String PANTHER_URL = "http://www.pantherdb.org/";
    
    public static final String MODEL_ELM_NAME = "model";
    public static Namespace SBML_NS = Namespace.getNamespace("http://www.sbml.org/sbml/level2");
    public static final Namespace SBML_NS_LEVEL1 = Namespace.getNamespace("http://www.sbml.org/sbml/level1");
    public static final Namespace CELL_DESIGNER_NS = Namespace.getNamespace("celldesigner",
            "http://www.sbml.org/2001/ns/celldesigner");
    public static final String ID_ATT_NAME = "id";
    public static final String PANTHER_SPACE_MASK = "_space_";
    public static final String EMPTY_SYNONYM_MARK = "Synonym not specified";
    public static final String EMPTY_LONG_NAME_MARK = "Long name not specified";
    // Seems a bug in the old version CellDesigner
    public static final String EMPTY_LONG_NAME_MARK_LEVEL1 = "Long Name:  Long name not specified";
    public static final String SYNONYM_LABEL = "Synonym";
    public static final String LONG_NAME_LABEL = "Long Name";
    public static final String ACCESSION_LABEL = "Accession";
    public static final String POSITION_TO_COMPARTMENT_ELM_NAME = "positionToCompartment";
    public static final String SPECIES_IDENTITY_ELM_NAME = "speciesIdentity";
    public static final String CLASS_ELM_NAME = "class";
    public static final String LIST_OF_SPECIES_ELM_NAME = "listOfSpecies";
    public static final String SPECIES_ELM_NAME = "species";
    public static final String SPECIES_ATT_NAME = "species";
    public static final String SPECIES_REFERENCE = "speciesReference";
    public static final String COMPARTMENT_ATT_NAME = "compartment";
    public static final String NAME_ATT_NAME = "name";
    public static final String NAME_ELM_NAME = "name";
    public static final String NOTES_ELM_NAME = "notes";
    public static final String HTML_ELM_NAME = "html";
    public static final String BODY_ELM_NAME = "body";
    public static final String ANNOTATION_ELM_NAME = "annotation";
    // These constants are for compartment definitions
    public static final String LIST_OF_COMPARTMENTS_ELM_NAME = "listOfCompartments";
    public static final String COMPARTMENT_ELM_NAME = "compartment";
    // These constants are for reaction definitions
    public static final String LIST_OF_REACTION_ELM_NAME = "listOfReactions";
    public static final String REACTION_ELM_NAME = "reaction";
    public static final String REVERSIBLE_ATT_NAME = "reversible";
    public static final String REACTION_TYPE_ELM_NAME = "reactionType";
    public static final String BASE_REACTANTS_ELM_NAME = "baseReactants";
    public static final String BASE_PRODUCTS_ELM_NAME = "baseProducts";
    public static final String LIST_OF_REACTANT_LINKS_ELM_NAME = "listOfReactantLinks";
    public static final String REACTANT_LINK_ELM_NAME = "reactantLink";
    public static final String REACTANT_ATT_NAME = "reactant";
    public static final String LIST_OF_PRODUCT_LINKS_ELM_NAME = "listOfProductLinks";
    public static final String PRODUCT_LINK_ELM_NAME = "productLink";
    public static final String PRODUCT_ATT_NAME = "product";
    public static final String LIST_OF_REACTANTS_ELM_NAME = "listOfReactants";
    public static final String LIST_OF_PRODUCTS_ELM_NAME = "listOfProducts";
    public static final String LIST_OF_MODIFIERS_ELM_NAME = "listOfModifiers";
    // For Reaction Modifications
    public static final String LIST_OF_MODIFICATION_ELM_NAME = "listOfModification";
    public static final String MODIFICATION_ELM_NAME = "modification";
    public static final String MODIFIFIER_ATT_NAME = "modifiers";
    public static final String modifierSpeciesReference = "modifierSpeciesReference";
    public static final String TYPE_ATT_NAME = "type";
    public static final String ALIASES_ATT_NAME = "aliases";
    public static final String ALIAS_ATT_NAME = "alias";
    // Types used for arrows that can be used for linking to both entities and reactions
    public static final String STATE_TRANSITION_TYPE = "STATE_TRANSITION";
    public static final String KNOWN_TRANSITION_OMITTED_TYPE = "KNOWN_TRANSITION_OMITTED";
    public static final String UNKNOWN_TRANSITION_TYPE = "UNKNOWN_TRANSITION";
    public static final String CATALYSIS_TYPE = "CATALYSIS";
    public static final String UNKNOWN_CATALYSIS_TYPE = "UNKNOWN_CATALYSIS";
    public static final String INHIBITION_TYPE = "INHIBITION";
    public static final String UNKNOWN_INHIBITION_TYPE = "UNKNOWN_INHIBITION";
    public static final String TRANSPORT_TYPE = "TRANSPORT";
    public static final String TRANSCRIPTIONAL_ACTIVATION_TYPE = "TRANSCRIPTIONAL_ACTIVATION";
    public static final String TRANSCRIPTIONAL_INHIBITION_TYPE = "TRANSCRIPTIONAL_INHIBITION";
    public static final String TRANSLATIONAL_ACTIVATION_TYPE = "TRANSLATIONAL_ACTIVATION";
    public static final String TRANSLATIONAL_INHIBITION_TYPE = "TRANSLATIONAL_INHIBITION";
    // These tags are for pathway notes
    // Namespace used for html. 
    public static final Namespace HTML_NS = Namespace.getNamespace("http://www.w3.org/1999/xhtml");
    public static final String WEB_SITE_LABEL = "Website";
    public static final String MEDLINE_LABEL = "Medline";
    public static final String PMID_LABEL = "PMID";
    public static final String FREE_TEXT_LABEL = "free text";
    public static final String EMPTY_PATHWAY_NOTE_MARK = "No pathway definition";
    public static final String PUBMED_LINE_MARK = "db=pubmed";
    public static final String PUBMED_ID_MARK = "list_uids=";
    // For compartments element
    public static final String DEFAULT_COMPARTMENT_NAME = "default";
    public static final String OUTSIDE_ATT_NAME = "outside";
    // For complex annotations
    public static final String HETERODIMER_IDENTITY_ELM_NAME = "heterodimerIdentity";
    public static final String LIST_OF_HETERODIMER_ENTRIES_ELM_NAME = "listOfHeterodimerEntries";
    public static final String HETERODIMER_ENTRY_ELM_NAME = "heterodimerEntry";
    public static final String PROTEIN_REFERENCE_ELM_NAME = "proteinReference";
    public static final String HOMODIMER_ELM_NAME = "homodimer";
    public static final String STATE_ELM_NAME = "state";
    // A list of type definitions
    public static final String PROTEIN_TYPE_NAME = "PROTEIN";
    public static final String GENE_TYPE_NAME = "GENE";
    public static final String RNA_TYPE_NAME = "RNA";
    public static final String ANTI_SENSE_RNA_TYPE_NAME = "ANTI_SENSE_RNA";
    public static final String PHENOTYPE_TYPE_NAME = "PHENOTYPE";
    public static final String ION_TYPE_NAME = "ION";
    public static final String SIMPLE_MOLECULE_TYPE_NAME = "SIMPLE_MOLECULE";
    public static final String UNKNOWN_TYPE_NAME = "UNKNOWN";
    public static final String DEGRADED_TYPE_NAME = "DEGRADED";
    public static final String GENERIC_TYPE_NAME = "GENERIC";
    // Customized, not used in CellDesigner
    public static final String COMPLEX_TYPE_NAME = "COMPLEX";
    // values of positionToCompartment
    public static final String INNER_SURFACE = "innerSurface";
    public static final String INSIDE = "inside";
    public static final String TRANSMEMBRANE = "transmembrane";
    public static final String OUTER_SURFACE = "outerSurface";
    // activity states
    public static final String ACTIVITY_ELM_NAME = "activity";
    public static final String INACTIVE_STATE = "inactive";
    public static final String ACTIVE_STATE = "active";
    // For Species modifications
    public static final String LIST_OF_MODIFICATION_RESIDUE_ELM_NAME = "listOfModificationResidues";
    public static final String MODIFICATION_RESIDUE_ELM_NAME = "modificationResidue";
    public static final String RESIDUE_ATT_NAME = "residue";
    public static final String STATE_ATT_NAME = "state";

    // for SBML2.2/2.3
    //public static final String LIST_OF_ALIASES_ELM_NAME = "listOfInnerAliases";
    public static final String LIST_OF_INNER_ALIASES_ELM_NAME = "listOfInnerAliases";
    public static final String INNERALIAS_ELM_NAME = "innerAlias";
    public static final String HETERODIMERENTRY_ATT_NAME = "heterodimerEntry";
    
    public static final String LIST_OF_INCLUDED_SPECIES_ELM_NAME = "listOfIncludedSpecies";
    public static final String COMPLEX_SPECIES_ELM_NAME = "complexSpecies";
    public static final String COMPARTMENT_ALIAS_ATT_NAME = "compartmentAlias";
    public static final String COMPLEX_SPECIES_ALIAS_ATT_NAME = "complexSpeciesAlias";
    
    // For model differentation and coordinates, height, width, x, y coordinates
    public static final String MODELVERSION_ELM_NAME = "modelVersion";
    public static final String BOUNDS_ELM_NAME = "bounds";
    public static final String USUALVIEW_ELM_NAME = "usualView";
    public static final String INNERPOSITION_ELM_NAME = "innerPosition";
    public static final String BOXSIZE_ELM_NAME = "boxSize";
    public static final String BOUNDS_TYPE_NAME = "bounds";
    public static final String H_ATT_NAME = "h";
    public static final String W_ATT_NAME = "w";
    public static final String X_ATT_NAME= "x";
    public static final String Y_ATT_NAME = "y";
    public static final String MODELDISPLAY_ELM_NAME = "modelDisplay";
    public static final String SIZEX_ATT_NAME = "sizeX";
    public static final String SIZEY_ATT_NAME = "sizeY";
    public static final String LIST_OF_COMPARTMENT_ALIASES_ELM_NAME = "listOfCompartmentAliases";
    public static final String ALIAS_ELM_NAME = "alias";
}
