/*
 * Created on Mar 20, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.List;

public class Interaction {
    private Long dbId;
    private int id;
    private Names names;
    private Xref xref;
    private List<Experiment> experimentList;
    private List<Participant> participantList;
    private List<InferredInteraction> inferredInteractionList;
    private OpenCV interactionType;
    private boolean modelled;
    private boolean intraMolecular;
    private boolean negative;
    private Availability availability;
    private List<ExperimentalConfidence> confidenceList;
    // No parameterList
    private List<Attribute> attributeList;
    private String imexId;
    
    public Interaction() {    
    }
    
    public String getImexId() {
        return imexId;
    }

    public void setImexId(String imexId) {
        this.imexId = imexId;
    }

    public void addExperiment(Experiment exp) {
        if (experimentList == null)
            experimentList = new ArrayList<Experiment>();
        experimentList.add(exp);
    }
    
    public void addParticipant(Participant participant) {
        if (participantList == null)
            participantList = new ArrayList<Participant>();
        participantList.add(participant);
    }
    
    public void addInferredInteraction(InferredInteraction interaction) {
        if (inferredInteractionList == null)
            inferredInteractionList = new ArrayList<InferredInteraction>();
        inferredInteractionList.add(interaction);
    }
    
    public void addConfidence(ExperimentalConfidence confidence) {
        if (confidenceList == null)
            confidenceList = new ArrayList<ExperimentalConfidence>();
        confidenceList.add(confidence);
    }
    
    public void addAttribute(Attribute attribute) {
        if (attributeList == null)
            attributeList = new ArrayList<Attribute>();
        attributeList.add(attribute);
    }
    
    public void setDbId(Long id) {
        this.dbId = id;
    }
    
    public Long getDbId() {
        return this.dbId;
    }
    
    public void setAvailability(Availability availability) {
        this.availability = availability;
    }
    
    public Availability getAvailability() {
        return this.availability;
    }

    public List<Attribute> getAttributeList() {
        return attributeList;
    }

    public void setAttributeList(List<Attribute> attributeList) {
        this.attributeList = attributeList;
    }

    public List<ExperimentalConfidence> getConfidenceList() {
        return confidenceList;
    }

    public void setConfidenceList(List<ExperimentalConfidence> confidenceList) {
        this.confidenceList = confidenceList;
    }

    public List<Experiment> getExperimentList() {
        return experimentList;
    }

    public void setExperimentList(List<Experiment> experimentList) {
        this.experimentList = experimentList;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public List<InferredInteraction> getInferredInteractionList() {
        return inferredInteractionList;
    }

    public void setInferredInteractionList(List<InferredInteraction> inferredInteractionList) {
        this.inferredInteractionList = inferredInteractionList;
    }

    public OpenCV getInteractionType() {
        return interactionType;
    }

    public void setInteractionType(OpenCV interactionType) {
        this.interactionType = interactionType;
    }

    public boolean getIntraMolecular() {
        return intraMolecular;
    }

    public void setIntraMolecular(boolean intraMolecular) {
        this.intraMolecular = intraMolecular;
    }

    public boolean getModelled() {
        return modelled;
    }

    public void setModelled(boolean modelled) {
        this.modelled = modelled;
    }

    public Names getNames() {
        return names;
    }

    public void setNames(Names names) {
        this.names = names;
    }

    public boolean getNegative() {
        return negative;
    }

    public void setNegative(boolean negative) {
        this.negative = negative;
    }

    public List<Participant> getParticipantList() {
        return participantList;
    }

    public void setParticipantList(List<Participant> participantList) {
        this.participantList = participantList;
    }

    public Xref getXref() {
        return xref;
    }

    public void setXref(Xref xref) {
        this.xref = xref;
    }
    
}
