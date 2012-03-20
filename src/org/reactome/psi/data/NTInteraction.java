/*
 * Created on Apr 27, 2006
 *
 */
package org.reactome.psi.data;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * This is a simplified PSI-MI interaction targetted to describe interactions among protein ids
 * with plain experiment names, interactor types.
 * @author guanming
 *
 */
public class NTInteraction {
    
    private Long dbId;
    // Use Set instead of List since position in the collection is not important
    // and save a little in hibernate mapping.
    private Set<NTInteractor> interactors;
    private Set<String> experimentTypes;
    private Set<String> interactionTypes;
    
    public NTInteraction() {
    }

    public Long getDbId() {
        return dbId;
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }

    public Set<String> getExperimentTypes() {
        return experimentTypes;
    }

    public void setExperimentTypes(Set<String> experimentTypes) {
        this.experimentTypes = experimentTypes;
    }
    
    public void addExperimentType(String type) {
        if (experimentTypes == null)
            experimentTypes = new HashSet<String>();
        experimentTypes.add(type);
    }

    public Set<NTInteractor> getInteractors() {
        return interactors;
    }

    public void setInteractors(Set<NTInteractor> interactors) {
        this.interactors = interactors;
    }
    
    public void addInteractor(NTInteractor interactor) {
        if (interactors == null)
            interactors = new HashSet<NTInteractor>();
        interactors.add(interactor);
    }
    
    public Set<String> getInteractionTypes() {
        return interactionTypes;
    }

    public void setInteractionType(Set<String> interactionTypes) {
        this.interactionTypes = interactionTypes;
    }
    
    public void addInteractionType(String type) {
        if (interactionTypes == null)
            interactionTypes = new HashSet<String>();
        interactionTypes.add(type);
    }

    public void exportToXML(StringBuilder xmlString) {
        xmlString.append("<interaction dbId=\"").append(dbId).append("\">");
        if (interactionTypes != null && interactionTypes.size() > 0) {
            xmlString.append("\n<interactionTypes>");
            for (Iterator it = interactionTypes.iterator(); it.hasNext();) {
                xmlString.append(it.next());
                if (it.hasNext())
                    xmlString.append(",");
            }
            xmlString.append("</interactionTypes>");
        }
        if (experimentTypes != null && experimentTypes.size() > 0) {
            xmlString.append("\n<experimentTypes>");
            for (Iterator<String> it = experimentTypes.iterator(); it.hasNext();) {
                xmlString.append(it.next());
                if (it.hasNext())
                    xmlString.append(",");
            }
            xmlString.append("</experimentTypes>");
        }
        if (interactors != null && interactors.size() > 0) {
            xmlString.append("\n<interactors>");
            for (NTInteractor interactor : interactors) {
                xmlString.append("\n");
                interactor.exportToXML(xmlString);
            }
            // A new line has been added at the end of the above loop.
            xmlString.append("\n</interactors>");
        }
        xmlString.append("\n</interaction>");
    }
}
