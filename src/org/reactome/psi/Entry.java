/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Entry {
    private Long dbId;
    private Source source;
    private List<Experiment> experimentList;
    private List<Interactor> interactorList;
    private List<Interaction> interactionList;
    private List<Attribute> attributeList;
    private List<Availability> availabilityList;
    
    public Entry() {
    }

    public void setDbId(Long dbId) {
        this.dbId = dbId;
    }
    
    public Long getDbId() {
        return this.dbId;
    }
    
    public List<Availability> getAvailabilityList() {
        return this.availabilityList;
    }
    
    public void setAvailabilityList(List<Availability> list) {
        this.availabilityList = list;
    }
    
    public void addAvailability(Availability availability) {
        if (availabilityList == null)
            availabilityList = new ArrayList<Availability>();
        availabilityList.add(availability);
    }
    
    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }
    
    public void addAttribute(Attribute att) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(att);
    }

    public List<Experiment> getExperimentList() {
        return experimentList;
    }

    public void setExperimentList(List<Experiment> experimentList) {
        this.experimentList = experimentList;
    }
    
    public void addExperiment(Experiment experiment) {
        if (experimentList == null)
            experimentList = new ArrayList<Experiment>();
        experimentList.add(experiment);
    }

    public List<Interaction> getInteractionList() {
        return interactionList;
    }

    public void setInteractionList(List<Interaction> interactionList) {
        this.interactionList = interactionList;
    }
    
    public void addInteraction(Interaction interaction) {
        if (interactionList == null)
            interactionList = new ArrayList<Interaction>();
        interactionList.add(interaction);
    }

    public List<Interactor> getInteractorList() {
        return interactorList;
    }

    public void setInteractorList(List<Interactor> interactorList) {
        this.interactorList = interactorList;
    }
    
    public void addInteractor(Interactor interactor) {
        if (interactorList == null)
            interactorList = new ArrayList<Interactor>();
        interactorList.add(interactor);
    }

    public Source getSource() {
        return source;
    }

    public void setSource(Source source) {
        this.source = source;
    }
}
