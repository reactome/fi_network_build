/*
 * Created on Mar 9, 2006
 *
 */
package org.reactome.panther;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.jdom.Element;

/**
 * This class is used to handle names for converting from Panther pathways to Reactome
 * pathways.
 * @author guanming
 */
public class PantherNameHandler {

    public PantherNameHandler() {
    }
    
    public String getSpeciesName(Element species) {
        String name = species.getAttributeValue(PantherConstants.NAME_ATT_NAME);
        return name;
    }
    
    public String getSpeciesID(Element species) {
        String id = species.getAttributeValue(PantherConstants.ID_ATT_NAME);
        if (id == null)
            id = getSpeciesName(species);
        return id;
    }
    
    public String cleanupComplexName(String name) {
   	 name = name.replace("Complex(", "");
   	 String pattern = "(_super_|_sub_)";
       name = name.replaceAll(pattern, "_");
       pattern = "(_endsuper_|_endsuper|_endsub_|_endsub)";
       name = name.replaceAll(pattern, "");
       name = name.replaceAll("(_br_|_space_)", " ");
       name = name.replaceAll("_minus_", "-");
       name = name.replaceAll("_plus_", "+");
       name = name.replaceAll("_slash_", "/");
       name = name.replaceAll("/", ":");
       // Use the encoded greek letters
       name = name.replaceAll("(_alpha_|_Alhpa_)", "_alpha");
       name = name.replaceAll("(_beta_|_Beta_)", "_beta");
       name = name.replaceAll("(_gamma_|_Gamma_)", "_gamma");
       name = name.replaceAll("(_delta_|_Delta_)", "_delta");
       name = name.trim();
       if (name.substring(name.length()-1).equals(")"))
      	 name = name.substring(0, name.length()-1);
       return name;
    }
    
    public String cleanupName(String name) {
        String pattern = "(_super_|_endsuper_|_endsuper|_sub_|_endsub_|_endsub)";
        name = name.replaceAll(pattern, "");
        name = name.replaceAll("(_br_|_space_)", " ");
        name = name.replaceAll("_minus_", "-");
        name = name.replaceAll("_plus_", "+");
        name = name.replaceAll("_slash_", "/");
        // Use the encoded greek letters
        name = name.replaceAll("(_alpha_|_Alhpa_|_alpha)", "alpha");
        name = name.replaceAll("(_beta_|_Beta_|_beta)", "beta");
        name = name.replaceAll("_Beta", "Beta");
        name = name.replaceAll("(_gamma_|_Gamma_)", "gamma");
        name = name.replaceAll("_delta_", "delta");
        return name.trim();
    }
    
    public List<String> getEntityNames(Element speciesElm, Map<String, Protein> proteinMap) {
        Set<String> nameSet = new HashSet<String>();
        // Names in level 2 pathways are useful
        // Check name attribute first. If there is no name,
        // try to get id
        String name = getSpeciesName(speciesElm);
        if (name == null)
            name = speciesElm.getAttributeValue(PantherConstants.ID_ATT_NAME);
        if (name != null)
            nameSet.add(name);
        // Figure out celldesigner:name, which is used for non-protein species.
        Element annotElm = speciesElm.getChild(PantherConstants.ANNOTATION_ELM_NAME,
                                               PantherConstants.SBML_NS);
        if (annotElm==null)
      	  annotElm = speciesElm.getChild(PantherConstants.ANNOTATION_ELM_NAME,
      			  									PantherConstants.CELL_DESIGNER_NS);
        Element speciesIdElm = annotElm.getChild(PantherConstants.SPECIES_IDENTITY_ELM_NAME,
                                                 PantherConstants.CELL_DESIGNER_NS);
        String proteinRef = null;
        String cellDesignerName = null;
        if (speciesIdElm != null) {
            cellDesignerName = speciesIdElm.getChildTextTrim(PantherConstants.NAME_ELM_NAME,
                                                        PantherConstants.CELL_DESIGNER_NS);
            if (cellDesignerName != null) {
                cellDesignerName = cleanupName(cellDesignerName);
                nameSet.add(cellDesignerName);
            }
            proteinRef = speciesIdElm.getChildTextTrim(PantherConstants.PROTEIN_REFERENCE_ELM_NAME,
                                                       PantherConstants.CELL_DESIGNER_NS);
        }
        if (proteinRef != null) {
            Protein protein = proteinMap.get(proteinRef);
            if (protein.getName() != null)
                nameSet.add(protein.getName());
        }
        List<String> nameList = new ArrayList<String>();
        if (nameSet.size() > 0) {
            nameList.addAll(nameSet);
            Collections.sort(nameList, new Comparator<String>() {
                public int compare(String name1, String name2) {
                    return name1.length() - name2.length();
                }
            });
            // Always place the cell designer name at the first
            if (cellDesignerName != null) {
                String firstName = nameList.get(0);
                if (!firstName.equals(cellDesignerName)) {
                    nameList.remove(cellDesignerName);
                    nameList.add(0, cellDesignerName);
                }
            }
        }
        return nameList;
    }
    
    @SuppressWarnings("unchecked")
    public void getNamesFromNotes(String[] lines, GKInstance entity) throws Exception {
        if (lines == null)
            return;
        List<String> names = entity.getAttributeValuesList(ReactomeJavaConstants.name);
        for (String line : lines) {
            line = line.trim();
            if (line.length() == 0)
                continue;
            if (line.startsWith(PantherConstants.LONG_NAME_LABEL)) {
                String longName = line.substring(PantherConstants.LONG_NAME_LABEL.length() + 1).trim();
                // Seems a bug in CellDesigner. See Long Name for calcium in Angiogenesis.xml.
                if (longName.startsWith(PantherConstants.LONG_NAME_LABEL)) {
                    longName = longName.substring(PantherConstants.LONG_NAME_LABEL.length() + 1).trim();
                }
                // Usually it is not the case
                if (!(longName.equals(PantherConstants.EMPTY_LONG_NAME_MARK) ||
                      longName.equals(PantherConstants.EMPTY_LONG_NAME_MARK_LEVEL1))) { // Name is set
                    // names might be null for small molecules in level 1 pathways
                    if (names == null || names.size() == 0)
                        entity.addAttributeValue(ReactomeJavaConstants.name, longName);
                    else if (!names.contains(longName)) {
                        String firstName = (String) names.get(0);
                        if (longName.length() < firstName.length())
                            names.add(0, longName); // Place long name at the first so that it can be picked
                                                    // for _displayName.
                        else
                            names.add(longName);
                    }
                }
            }
            else if (line.startsWith(PantherConstants.SYNONYM_LABEL)) {
                String synonym = line.substring(PantherConstants.SYNONYM_LABEL.length() + 1).trim();
                if (!synonym.equals(PantherConstants.EMPTY_SYNONYM_MARK)) {
                    // Delimited by ","
                    String[] nameArray = synonym.split(",");
                    for (String name : nameArray) {
                        // Synonym might be the same as long name: CREB in Wnt Signaling Pathway.xml.
                        if (names != null && !names.contains(name))
                            entity.addAttributeValue(ReactomeJavaConstants.name,
                                    name.trim());
                    }
                }
            }
        }
    }
    
    public void assignDegradedEntityNames(List<GKInstance> degradedEntities, 
                                          Map<String, GKInstance> idToInstanceMap) throws Exception {
        if (degradedEntities == null || degradedEntities.size() == 0)
            return;
        for (GKInstance entity : degradedEntities) {
            String name = (String) entity.getAttributeValue(ReactomeJavaConstants.name);
            int index = name.indexOf("_");
            if (index < 0)
                continue;
            String speciesId = name.substring(0, index);
            GKInstance input = idToInstanceMap.get(speciesId);
            if (input != null) {
                String inputName = (String) input.getAttributeValue(ReactomeJavaConstants.name);
                String newName = inputName + name.substring(index);
                entity.setAttributeValueNoCheck(ReactomeJavaConstants.name, newName);
                InstanceDisplayNameGenerator.setDisplayName(entity);
            }
        }
    }
    
    public void validateComplexNames(Map<String, GKInstance> idToInstanceMap) throws Exception {
        Collection<GKInstance> instances = idToInstanceMap.values();
        for (GKInstance instance : instances) {
            if (instance.getSchemClass().isa(ReactomeJavaConstants.Complex))
                validateComplexName(instance);
        }
    }
    
    private void validateComplexName(GKInstance complex) throws Exception {
        // Create a name dynamically if no name defined for complex
        String name = (String) complex.getAttributeValue(ReactomeJavaConstants.name);
        if (name == null) {
            List comps = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
            if (comps != null && comps.size() > 0) {
                StringBuilder builder = new StringBuilder();
                builder.append("Complex(");
                for (Iterator it = comps.iterator(); it.hasNext();) {
                    GKInstance comp = (GKInstance) it.next();
                    String compName = (String) comp.getAttributeValue(ReactomeJavaConstants.name);
                    if (compName != null)
                        builder.append(compName);
                    if (it.hasNext())
                        builder.append("/");
                }
                builder.append(")");
                complex.addAttributeValue(ReactomeJavaConstants.name, builder.toString());
                InstanceDisplayNameGenerator.setDisplayName(complex);
            }
        }
    }
}
